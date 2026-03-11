! libadobmd_parser.f90 - Fixed memory management
! Compile: gfortran -shared -fPIC -O3 -o libadobmd_parser.so libadobmd_parser.f90

module adobmd_parser
    use iso_c_binding
    implicit none
    
    ! Constants
    integer, parameter :: MAX_LINE_LEN = 256
    integer, parameter :: MAX_FILENAME_LEN = 256
    
    ! Atom structure (C-compatible)
    type, bind(C) :: atom_t
        integer(c_int) :: id
        integer(c_int) :: type_id
        integer(c_int) :: molecule
        character(c_char) :: element(3)
        real(c_double) :: x, y, z
        real(c_double) :: charge
        integer(c_int) :: is_qm
    end type atom_t
    
    ! Bond structure
    type, bind(C) :: bond_t
        integer(c_int) :: id
        integer(c_int) :: type_id
        integer(c_int) :: atom1
        integer(c_int) :: atom2
        integer(c_int) :: order
    end type bond_t
    
    ! Parser result
    type, bind(C) :: parse_result_t
        integer(c_int) :: natoms
        integer(c_int) :: nbonds
        integer(c_int) :: nqm
        real(c_double) :: box_lo(3)
        real(c_double) :: box_hi(3)
        integer(c_int) :: has_box
        type(c_ptr) :: atoms
        type(c_ptr) :: bonds
        character(c_char) :: filename(MAX_FILENAME_LEN)
    end type parse_result_t
    
contains
    
    ! Helper: Convert C string to Fortran string
    subroutine c_f_string(c_str, f_str)
        character(c_char), intent(in) :: c_str(*)
        character(len=*), intent(out) :: f_str
        integer :: i
        
        f_str = ''
        i = 1
        do while (c_str(i) /= c_null_char .and. i < len(f_str))
            f_str(i:i) = c_str(i)
            i = i + 1
        end do
    end subroutine c_f_string
    
    ! Helper: Convert Fortran string to C string
    subroutine f_c_string(f_str, c_str)
        character(len=*), intent(in) :: f_str
        character(c_char), intent(out) :: c_str(*)
        integer :: i, len_str
        
        len_str = len_trim(f_str)
        do i = 1, len_str
            c_str(i) = f_str(i:i)
        end do
        c_str(len_str + 1) = c_null_char
    end subroutine f_c_string
    
    ! Parse header
    subroutine parse_header(unit, natoms, nbonds, natom_types, has_box, box_lo, box_hi)
        integer, intent(in) :: unit
        integer, intent(out) :: natoms, nbonds, natom_types
        logical, intent(out) :: has_box
        real(8), intent(out) :: box_lo(3), box_hi(3)
        
        character(len=MAX_LINE_LEN) :: line
        integer :: ios
        logical :: reading_header
        
        natoms = 0; nbonds = 0; natom_types = 0
        has_box = .false.
        box_lo = 0.0d0; box_hi = 0.0d0
        reading_header = .true.
        
        do while (reading_header)
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            line = adjustl(line)
            if (line(1:1) == '#') cycle
            if (len_trim(line) == 0) cycle
            
            if (index(line, 'atoms') > 0) then
                read(line, *) natoms
            else if (index(line, 'bonds') > 0) then
                read(line, *) nbonds
            else if (index(line, 'atom types') > 0) then
                read(line, *) natom_types
            else if (index(line, 'xlo xhi') > 0) then
                read(line, *) box_lo(1), box_hi(1)
                has_box = .true.
            else if (index(line, 'ylo yhi') > 0) then
                read(line, *) box_lo(2), box_hi(2)
            else if (index(line, 'zlo zhi') > 0) then
                read(line, *) box_lo(3), box_hi(3)
            else if (line == 'Masses') then
                reading_header = .false.
                exit
            end if
        end do
        
        if (reading_header) then
            ! Rewind if we didn't find Masses
            backspace(unit)
        end if
    end subroutine parse_header
    
    ! Parse type_elements from Masses section
    subroutine parse_type_elements(unit, natom_types, type_elements)
        integer, intent(in) :: unit, natom_types
        character(len=2), intent(out) :: type_elements(:)
        
        character(len=MAX_LINE_LEN) :: line
        integer :: i, type_id, ios, idx
        real(8) :: mass
        integer :: comment_pos
        
        ! Skip to Masses section
        rewind(unit)
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (line == 'Masses') exit
        end do
        
        ! Skip blank lines
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
        end do
        backspace(unit)
        
        ! Read masses
        do i = 1, natom_types
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            ! Extract type_id and mass
            read(line, *) type_id, mass
            
            ! Extract element from comment
            type_elements(type_id) = 'X'
            comment_pos = index(line, '#')
            if (comment_pos > 0) then
                do idx = comment_pos + 1, len_trim(line)
                    if (line(idx:idx) >= 'A' .and. line(idx:idx) <= 'Z') then
                        if (idx + 1 <= len_trim(line) .and. &
                            line(idx+1:idx+1) >= 'a' .and. line(idx+1:idx+1) <= 'z') then
                            type_elements(type_id) = line(idx:idx+1)
                        else
                            type_elements(type_id) = line(idx:idx)
                        end if
                        exit
                    end if
                end do
            end if
        end do
    end subroutine parse_type_elements
    
    ! Parse atoms section
    subroutine parse_atoms(unit, natoms, atoms, type_elements, qm_indices, nqm)
        integer, intent(in) :: unit, natoms
        type(atom_t), intent(out) :: atoms(:)
        character(len=2), intent(in) :: type_elements(:)
        integer, intent(out) :: qm_indices(:)
        integer, intent(out) :: nqm
        
        character(len=MAX_LINE_LEN) :: line
        integer :: i, ios, atom_id, type_id, molecule
        real(8) :: x, y, z, charge
        character(len=2) :: element
        logical :: is_qm
        
        ! Skip to Atoms section
        rewind(unit)
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (line == 'Atoms') exit
        end do
        
        ! Skip blank lines
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
        end do
        backspace(unit)
        
        ! Read atoms
        nqm = 0
        do i = 1, natoms
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            ! Parse atom line
            read(line, *) atom_id, molecule, type_id, charge, x, y, z
            
            ! Check for region tag
            is_qm = .false.
            if (index(line, 'QM') > 0 .or. index(line, 'Q') > 0) then
                is_qm = .true.
                nqm = nqm + 1
                qm_indices(nqm) = atom_id
            end if
            
            ! Get element from type
            element = type_elements(type_id)
            
            ! Fill atom structure
            atoms(i)%id = atom_id
            atoms(i)%type_id = type_id
            atoms(i)%molecule = molecule
            atoms(i)%element = transfer(trim(element)//c_null_char, atoms(i)%element)
            atoms(i)%x = x
            atoms(i)%y = y
            atoms(i)%z = z
            atoms(i)%charge = charge
            atoms(i)%is_qm = merge(1, 0, is_qm)
        end do
    end subroutine parse_atoms
    
    ! Parse bonds section
    subroutine parse_bonds(unit, nbonds, bonds)
        integer, intent(in) :: unit, nbonds
        type(bond_t), intent(out) :: bonds(:)
        
        character(len=MAX_LINE_LEN) :: line
        integer :: i, ios, bond_id, bond_type, atom1, atom2
        
        if (nbonds == 0) return
        
        ! Skip to Bonds section
        rewind(unit)
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (line == 'Bonds') exit
        end do
        
        ! Skip blank lines
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
        end do
        backspace(unit)
        
        ! Read bonds
        do i = 1, nbonds
            read(unit, *, iostat=ios) bond_id, bond_type, atom1, atom2
            if (ios /= 0) exit
            
            bonds(i)%id = bond_id
            bonds(i)%type_id = bond_type
            bonds(i)%atom1 = atom1
            bonds(i)%atom2 = atom2
            bonds(i)%order = bond_type
        end do
    end subroutine parse_bonds
    
    ! Main parsing function
    subroutine parse_file(filename_c, result_c) bind(C, name="parse_file")
        character(c_char), intent(in) :: filename_c(*)
        type(parse_result_t), intent(out) :: result_c
        
        character(len=256) :: filename
        integer :: unit, natoms, nbonds, natom_types, nqm, ios
        real(8) :: box_lo(3), box_hi(3)
        logical :: has_box
        type(atom_t), pointer :: atoms(:)
        type(bond_t), pointer :: bonds(:)
        character(len=2), allocatable :: type_elements(:)
        integer, allocatable :: qm_indices(:)
        
        ! Initialize result
        result_c%natoms = 0
        result_c%nbonds = 0
        result_c%atoms = c_null_ptr
        result_c%bonds = c_null_ptr
        
        ! Convert C string to Fortran string
        call c_f_string(filename_c, filename)
        
        ! Store filename
        call f_c_string(trim(filename), result_c%filename)
        
        ! Open file
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            result_c%natoms = -1
            return
        end if
        
        ! Parse header
        call parse_header(unit, natoms, nbonds, natom_types, has_box, box_lo, box_hi)
        
        ! Allocate memory
        allocate(atoms(natoms))
        allocate(bonds(nbonds))
        allocate(type_elements(natom_types))
        allocate(qm_indices(natoms))
        
        ! Parse type elements
        call parse_type_elements(unit, natom_types, type_elements)
        
        ! Parse atoms
        call parse_atoms(unit, natoms, atoms, type_elements, qm_indices, nqm)
        
        ! Parse bonds
        if (nbonds > 0) then
            call parse_bonds(unit, nbonds, bonds)
        end if
        
        close(unit)
        
        ! Fill result structure
        result_c%natoms = natoms
        result_c%nbonds = nbonds
        result_c%nqm = nqm
        result_c%box_lo = box_lo
        result_c%box_hi = box_hi
        result_c%has_box = merge(1, 0, has_box)
        result_c%atoms = c_loc(atoms)
        result_c%bonds = c_loc(bonds)
        
        ! Deallocate temporary arrays (atoms and bonds are kept, they're now managed by C)
        deallocate(type_elements)
        deallocate(qm_indices)
        
    end subroutine parse_file
    
    ! Free memory function - CRITICAL for avoiding double free
    subroutine free_result(result_c) bind(C, name="free_result")
        type(parse_result_t), intent(inout) :: result_c
        type(atom_t), pointer :: atoms(:)
        type(bond_t), pointer :: bonds(:)
        
        ! Only deallocate if pointers are associated
        if (c_associated(result_c%atoms)) then
            call c_f_pointer(result_c%atoms, atoms, [result_c%natoms])
            deallocate(atoms)
            result_c%atoms = c_null_ptr
        end if
        
        if (c_associated(result_c%bonds)) then
            call c_f_pointer(result_c%bonds, bonds, [result_c%nbonds])
            deallocate(bonds)
            result_c%bonds = c_null_ptr
        end if
        
        ! Reset counts
        result_c%natoms = 0
        result_c%nbonds = 0
        result_c%nqm = 0
        
    end subroutine free_result
    
    ! Get version information
    subroutine get_version(version_str) bind(C, name="get_version")
        character(c_char), intent(out) :: version_str(*)
        call f_c_string('ADOBMD Fortran Parser v1.0.1 (fixed memory)', version_str)
    end subroutine get_version

end module adobmd_parser