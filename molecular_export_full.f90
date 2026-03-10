! molecular_export_full.f90 - Complete molecular export with Schlegel projections
! Exports: XYZ, SDF, PDB, MOL2, CIF, GRO + XY/XZ/YZ plane projections
! Compile: gfortran -shared -fPIC -O3 -o libmolexport.so molecular_export_full.f90

module molecular_export
    use iso_c_binding
    implicit none
    
    ! Constants
    integer, parameter :: MAX_LINE_LEN = 256
    integer, parameter :: MAX_FILENAME_LEN = 256
    integer, parameter :: MAX_ELEMENT_LEN = 3
    real(8), parameter :: PI = 3.14159265358979323846d0
    real(8), parameter :: DEG_TO_RAD = PI / 180.0d0
    real(8), parameter :: RAD_TO_DEG = 180.0d0 / PI
    
    ! Atom structure
    type, bind(C) :: export_atom_t
        integer(c_int) :: id
        integer(c_int) :: type_id
        integer(c_int) :: molecule
        character(c_char) :: element(3)
        real(c_double) :: x, y, z
        real(c_double) :: charge
        integer(c_int) :: is_qm
        ! Projections for Schlegel diagrams
        real(c_double) :: xy_x, xy_y      ! XY plane (top view)
        real(c_double) :: xz_x, xz_z      ! XZ plane (side view)
        real(c_double) :: yz_y, yz_z      ! YZ plane (side view)
        ! Scaled projections for visualization
        real(c_double) :: schlegel_x, schlegel_y  ! Main Schlegel view
    end type export_atom_t
    
    ! Bond structure
    type, bind(C) :: export_bond_t
        integer(c_int) :: id
        integer(c_int) :: type_id
        integer(c_int) :: atom1
        integer(c_int) :: atom2
        integer(c_int) :: order
    end type export_bond_t
    
    ! Export result
    type, bind(C) :: export_result_t
        integer(c_int) :: natoms
        integer(c_int) :: nbonds
        integer(c_int) :: nqm
        real(c_double) :: box_lo(3)
        real(c_double) :: box_hi(3)
        integer(c_int) :: has_box
        type(c_ptr) :: atoms
        type(c_ptr) :: bonds
        character(c_char) :: filename(MAX_FILENAME_LEN)
        ! Schlegel parameters
        real(c_double) :: projection_point(3)
        real(c_double) :: view_angle
        real(c_double) :: scale_factor
    end type export_result_t
    
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
    
    ! Calculate rotation matrix to align vector with z-axis
    subroutine rotation_matrix_to_z_axis(vec, rot_mat)
        real(8), intent(in) :: vec(3)
        real(8), intent(out) :: rot_mat(3,3)
        
        real(8) :: v_norm, v1(3), theta, phi, cos_theta, sin_theta, cos_phi, sin_phi
        
        v_norm = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
        if (v_norm < 1.0d-10) then
            rot_mat = reshape([1.0d0, 0.0d0, 0.0d0, &
                               0.0d0, 1.0d0, 0.0d0, &
                               0.0d0, 0.0d0, 1.0d0], [3,3])
            return
        end if
        
        ! Normalize vector
        v1 = vec / v_norm
        
        ! Calculate spherical angles
        theta = acos(v1(3))
        phi = atan2(v1(2), v1(1))
        
        cos_theta = cos(theta)
        sin_theta = sin(theta)
        cos_phi = cos(phi)
        sin_phi = sin(phi)
        
        ! Rotation matrix: Rz(phi) * Ry(theta)
        rot_mat(1,1) = cos_theta * cos_phi
        rot_mat(1,2) = -sin_phi
        rot_mat(1,3) = sin_theta * cos_phi
        
        rot_mat(2,1) = cos_theta * sin_phi
        rot_mat(2,2) = cos_phi
        rot_mat(2,3) = sin_theta * sin_phi
        
        rot_mat(3,1) = -sin_theta
        rot_mat(3,2) = 0.0d0
        rot_mat(3,3) = cos_theta
    end subroutine rotation_matrix_to_z_axis
    
    ! Apply rotation to a vector
    subroutine apply_rotation(rot_mat, vec, vec_out)
        real(8), intent(in) :: rot_mat(3,3), vec(3)
        real(8), intent(out) :: vec_out(3)
        
        vec_out(1) = rot_mat(1,1)*vec(1) + rot_mat(1,2)*vec(2) + rot_mat(1,3)*vec(3)
        vec_out(2) = rot_mat(2,1)*vec(1) + rot_mat(2,2)*vec(2) + rot_mat(2,3)*vec(3)
        vec_out(3) = rot_mat(3,1)*vec(1) + rot_mat(3,2)*vec(2) + rot_mat(3,3)*vec(3)
    end subroutine apply_rotation
    
    ! Calculate perspective projection for Schlegel diagram
    subroutine perspective_projection(point, projection_point, scale, proj_point)
        real(8), intent(in) :: point(3)
        real(8), intent(in) :: projection_point(3)  ! Point from which we project
        real(8), intent(in) :: scale                 ! Scaling factor
        real(8), intent(out) :: proj_point(2)
        
        real(8) :: z_dist, factor
        
        z_dist = projection_point(3) - point(3)
        if (abs(z_dist) < 1.0d-10) then
            proj_point(1) = point(1) * scale
            proj_point(2) = point(2) * scale
        else
            factor = (projection_point(3) - minval([point(3)])) / z_dist * scale
            proj_point(1) = point(1) * factor
            proj_point(2) = point(2) * factor
        end if
    end subroutine perspective_projection
    
    ! Calculate cone projection for Schlegel diagram
    subroutine cone_projection(point, projection_point, angle, proj_point)
        real(8), intent(in) :: point(3)
        real(8), intent(in) :: projection_point(3)
        real(8), intent(in) :: angle  ! Cone angle in degrees
        real(8), intent(out) :: proj_point(2)
        
        real(8) :: r, z_dist, factor, angle_rad
        
        angle_rad = angle * DEG_TO_RAD
        z_dist = projection_point(3) - point(3)
        r = sqrt(point(1)**2 + point(2)**2)
        
        if (r < 1.0d-10) then
            proj_point = 0.0d0
        else
            factor = tan(angle_rad) * z_dist / r
            proj_point(1) = point(1) * factor
            proj_point(2) = point(2) * factor
        end if
    end subroutine cone_projection
    
    ! Calculate all Schlegel projections for atoms
    subroutine calculate_schlegel_projections(natoms, atoms, projection_point, view_angle)
        integer, intent(in) :: natoms
        type(export_atom_t), intent(inout) :: atoms(:)
        real(8), intent(in) :: projection_point(3)
        real(8), intent(in) :: view_angle
        
        integer :: i
        real(8) :: point(3), proj_xy(2), proj_xz(2), proj_yz(2)
        
        do i = 1, natoms
            point = [atoms(i)%x, atoms(i)%y, atoms(i)%z]
            
            ! Basic projections (orthographic)
            atoms(i)%xy_x = point(1)
            atoms(i)%xy_y = point(2)
            
            atoms(i)%xz_x = point(1)
            atoms(i)%xz_z = point(3)
            
            atoms(i)%yz_y = point(2)
            atoms(i)%yz_z = point(3)
            
            ! Perspective projection for main Schlegel view
            call perspective_projection(point, projection_point, 1.0d0, proj_xy)
            atoms(i)%schlegel_x = proj_xy(1)
            atoms(i)%schlegel_y = proj_xy(2)
        end do
    end subroutine calculate_schlegel_projections
    
    ! Export XYZ file
    subroutine export_xyz(filename, natoms, atoms, with_qm_marker)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: natoms
        type(export_atom_t), intent(in) :: atoms(:)  ! Changed from assumed-size
        logical, optional :: with_qm_marker
        
        integer :: unit, i
        logical :: qm_marker
        character(len=10) :: element_str
        
        qm_marker = .false.
        if (present(with_qm_marker)) qm_marker = with_qm_marker
        
        open(newunit=unit, file=filename, status='replace', action='write')
        
        write(unit, '(I0)') natoms
        if (qm_marker) then
            write(unit, '(A)') 'ADOBMD Export with QM/MM regions'
        else
            write(unit, '(A)') 'ADOBMD Export'
        end if
        
        do i = 1, natoms
            call c_f_string(atoms(i)%element, element_str)
            write(unit, '(A2,3F15.8)') trim(element_str), atoms(i)%x, atoms(i)%y, atoms(i)%z
        end do
        
        close(unit)
    end subroutine export_xyz
    
    ! Export PDB file
    subroutine export_pdb(filename, natoms, atoms, bonds, nbonds)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: natoms, nbonds
        type(export_atom_t), intent(in) :: atoms(:)  ! Changed from assumed-size
        type(export_bond_t), intent(in) :: bonds(:)  ! Changed from assumed-size
        
        integer :: unit, i
        character(len=3) :: res_name
        character(len=2) :: element_str
        
        open(newunit=unit, file=filename, status='replace', action='write')
        
        ! Header
        write(unit, '(A)') 'HEADER    ADOBMD MOLECULAR STRUCTURE'
        write(unit, '(A,I6,A,I6)') 'REMARK    ATOMS: ', natoms, '  BONDS: ', nbonds
        write(unit, '(A)') 'REMARK    Generated by ADOBMD Molecular Exporter'
        
        ! Atoms
        do i = 1, natoms
            call c_f_string(atoms(i)%element, element_str)
            res_name = 'MOL'
            
            ! PDB format: fixed column positions
            write(unit, '(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2)') &
                'ATOM  ', i, adjustr(element_str)//'  ', ' ', &
                res_name, 'A', atoms(i)%molecule, ' ', &
                atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                1.0d0, 0.0d0, adjustl(element_str)
        end do
        
        ! Bonds
        if (nbonds > 0) then
            do i = 1, nbonds
                write(unit, '(A6,4I5)') 'CONECT', bonds(i)%atom1, bonds(i)%atom2, &
                    bonds(i)%atom1, bonds(i)%atom2
            end do
        end if
        
        write(unit, '(A)') 'END'
        close(unit)
    end subroutine export_pdb
    
    ! Export SDF file
    subroutine export_sdf(filename, natoms, nbonds, atoms, bonds)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: natoms, nbonds
        type(export_atom_t), intent(in) :: atoms(:)  ! Changed from assumed-size
        type(export_bond_t), intent(in) :: bonds(:)  ! Changed from assumed-size
        
        integer :: unit, i
        character(len=2) :: element_str
        character(len=20) :: atom_line
        
        open(newunit=unit, file=filename, status='replace', action='write')
        
        ! Header
        write(unit, '(A)') 'ADOBMD'
        write(unit, '(A)') '  ADOBMD Export'
        write(unit, '(A)') ''
        
        ! Counts line - fixed format
        write(unit, '(I3,I3,6I3)') natoms, nbonds, 0, 0, 0, 0, 0, 0
        
        ! Atoms - fixed format string
        do i = 1, natoms
            call c_f_string(atoms(i)%element, element_str)
            ! Format: 3F10.4, then element, then 3I2
            write(unit, '(3F10.4,3X,A2,3I3)') &
                atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                element_str, 0, 0, 0
        end do
        
        ! Bonds
        do i = 1, nbonds
            write(unit, '(3I3,3I3)') bonds(i)%atom1, bonds(i)%atom2, &
                bonds(i)%order, 0, 0, 0
        end do
        
        write(unit, '(A)') 'M  END'
        close(unit)
    end subroutine export_sdf
    
    ! Export MOL2 file
    subroutine export_mol2(filename, natoms, nbonds, atoms, bonds)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: natoms, nbonds
        type(export_atom_t), intent(in) :: atoms(:)  ! Changed from assumed-size
        type(export_bond_t), intent(in) :: bonds(:)  ! Changed from assumed-size
        
        integer :: unit, i
        character(len=2) :: element_str
        character(len=5) :: atom_type
        
        open(newunit=unit, file=filename, status='replace', action='write')
        
        ! Molecule section
        write(unit, '(A)') '@<TRIPOS>MOLECULE'
        write(unit, '(A)') 'ADOBMD'
        write(unit, '(I6,I6)') natoms, nbonds
        write(unit, '(A)') 'SMALL'
        write(unit, '(A)') 'NO_CHARGES'
        write(unit, '(A)') ''
        
        ! Atom section
        write(unit, '(A)') '@<TRIPOS>ATOM'
        do i = 1, natoms
            call c_f_string(atoms(i)%element, element_str)
            atom_type = element_str
            ! Fixed format
            write(unit, '(I6,1X,A4,3F12.6,1X,A5,1X,I2,1X,A4,F12.6)') &
                i, element_str, atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                atom_type, 1, 'MOL', atoms(i)%charge
        end do
        
        ! Bond section
        if (nbonds > 0) then
            write(unit, '(A)') '@<TRIPOS>BOND'
            do i = 1, nbonds
                write(unit, '(I6,I6,I6,I6)') i, bonds(i)%atom1, bonds(i)%atom2, bonds(i)%order
            end do
        end if
        
        close(unit)
    end subroutine export_mol2
    
    ! Export CIF file
    subroutine export_cif(filename, natoms, atoms, has_box, box_lo, box_hi)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: natoms
        type(export_atom_t), intent(in) :: atoms(:)  ! Changed from assumed-size
        logical, intent(in) :: has_box
        real(8), intent(in) :: box_lo(3), box_hi(3)
        
        integer :: unit, i
        character(len=2) :: element_str
        real(8) :: a, b, c
        character(len=10) :: label
        
        open(newunit=unit, file=filename, status='replace', action='write')
        
        write(unit, '(A)') 'data_adobmd'
        write(unit, '(A)') '_audit_creation_method ADOBMD Export'
        write(unit, '(A)') ''
        
        if (has_box) then
            a = box_hi(1) - box_lo(1)
            b = box_hi(2) - box_lo(2)
            c = box_hi(3) - box_lo(3)
            write(unit, '(A,F10.4)') '_cell_length_a ', a
            write(unit, '(A,F10.4)') '_cell_length_b ', b
            write(unit, '(A,F10.4)') '_cell_length_c ', c
            write(unit, '(A)') '_cell_angle_alpha 90.0'
            write(unit, '(A)') '_cell_angle_beta 90.0'
            write(unit, '(A)') '_cell_angle_gamma 90.0'
            write(unit, '(A)') ''
        end if
        
        write(unit, '(A)') 'loop_'
        write(unit, '(A)') '_atom_site_label'
        write(unit, '(A)') '_atom_site_type_symbol'
        write(unit, '(A)') '_atom_site_fract_x'
        write(unit, '(A)') '_atom_site_fract_y'
        write(unit, '(A)') '_atom_site_fract_z'
        write(unit, '(A)') '_atom_site_charge'
        
        do i = 1, natoms
            call c_f_string(atoms(i)%element, element_str)
            write(label, '(A,I0)') 'A', i
            write(unit, '(A6,1X,A2,3F12.6,F8.3)') &
                label, element_str, &
                atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                atoms(i)%charge
        end do
        
        close(unit)
    end subroutine export_cif
    
    ! Export GRO file
    subroutine export_gro(filename, natoms, atoms, has_box, box_lo, box_hi)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: natoms
        type(export_atom_t), intent(in) :: atoms(:)  ! Changed from assumed-size
        logical, intent(in) :: has_box
        real(8), intent(in) :: box_lo(3), box_hi(3)
        
        integer :: unit, i
        character(len=2) :: element_str
        real(8) :: a, b, c
        character(len=5) :: res_name
        
        open(newunit=unit, file=filename, status='replace', action='write')
        
        write(unit, '(A)') 'ADOBMD Export'
        write(unit, '(I5)') natoms
        write(unit, '(A)') ''
        
        res_name = 'MOL'
        do i = 1, natoms
            call c_f_string(atoms(i)%element, element_str)
            ! GRO format: residue number(5), residue name(5), atom name(5), atom number(5), positions(3F8.3)
            write(unit, '(I5,2X,A5,2X,A5,2X,I5,3F8.3)') &
                1, res_name, element_str, i, &
                atoms(i)%x/10.0d0, atoms(i)%y/10.0d0, atoms(i)%z/10.0d0
        end do
        
        if (has_box) then
            a = (box_hi(1) - box_lo(1)) / 10.0d0
            b = (box_hi(2) - box_lo(2)) / 10.0d0
            c = (box_hi(3) - box_lo(3)) / 10.0d0
            write(unit, '(3F10.5)') a, b, c
        end if
        
        close(unit)
    end subroutine export_gro
    
    ! Export Schlegel projections to separate files
    subroutine export_schlegel_projections(base_filename, natoms, atoms, projection_point, view_angle)
        character(len=*), intent(in) :: base_filename
        integer, intent(in) :: natoms
        type(export_atom_t), intent(in) :: atoms(:)  ! Changed from assumed-size
        real(8), intent(in) :: projection_point(3)
        real(8), intent(in) :: view_angle
        
        integer :: unit, i
        character(len=MAX_FILENAME_LEN) :: filename
        character(len=2) :: element_str
        
        ! Export XY plane projection
        filename = trim(base_filename) // '_xy.txt'
        open(newunit=unit, file=filename, status='replace', action='write')
        write(unit, '(A)') '# XY Plane Projection (top view)'
        write(unit, '(A,3F10.4)') '# Projection point: ', projection_point
        write(unit, '(A,F10.4)') '# View angle: ', view_angle
        write(unit, '(A)') '# ID  Element  X  Y  Z  Proj_X  Proj_Y'
        do i = 1, natoms
            call c_f_string(atoms(i)%element, element_str)
            write(unit, '(I5,1X,A2,5F10.4)') &
                atoms(i)%id, element_str, &
                atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                atoms(i)%xy_x, atoms(i)%xy_y
        end do
        close(unit)
        
        ! Export XZ plane projection
        filename = trim(base_filename) // '_xz.txt'
        open(newunit=unit, file=filename, status='replace', action='write')
        write(unit, '(A)') '# XZ Plane Projection (side view)'
        write(unit, '(A,3F10.4)') '# Projection point: ', projection_point
        write(unit, '(A,F10.4)') '# View angle: ', view_angle
        write(unit, '(A)') '# ID  Element  X  Y  Z  Proj_X  Proj_Z'
        do i = 1, natoms
            call c_f_string(atoms(i)%element, element_str)
            write(unit, '(I5,1X,A2,5F10.4)') &
                atoms(i)%id, element_str, &
                atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                atoms(i)%xz_x, atoms(i)%xz_z
        end do
        close(unit)
        
        ! Export YZ plane projection
        filename = trim(base_filename) // '_yz.txt'
        open(newunit=unit, file=filename, status='replace', action='write')
        write(unit, '(A)') '# YZ Plane Projection (side view)'
        write(unit, '(A,3F10.4)') '# Projection point: ', projection_point
        write(unit, '(A,F10.4)') '# View angle: ', view_angle
        write(unit, '(A)') '# ID  Element  X  Y  Z  Proj_Y  Proj_Z'
        do i = 1, natoms
            call c_f_string(atoms(i)%element, element_str)
            write(unit, '(I5,1X,A2,5F10.4)') &
                atoms(i)%id, element_str, &
                atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                atoms(i)%yz_y, atoms(i)%yz_z
        end do
        close(unit)
        
        ! Export main Schlegel view (perspective projection)
        filename = trim(base_filename) // '_schlegel.txt'
        open(newunit=unit, file=filename, status='replace', action='write')
        write(unit, '(A)') '# Schlegel Diagram (perspective projection)'
        write(unit, '(A,3F10.4)') '# Projection point: ', projection_point
        write(unit, '(A,F10.4)') '# View angle: ', view_angle
        write(unit, '(A)') '# ID  Element  X  Y  Z  Schlegel_X  Schlegel_Y'
        do i = 1, natoms
            call c_f_string(atoms(i)%element, element_str)
            write(unit, '(I5,1X,A2,5F10.4)') &
                atoms(i)%id, element_str, &
                atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                atoms(i)%schlegel_x, atoms(i)%schlegel_y
        end do
        close(unit)
        
        ! Export QM/MM colored version if applicable
        if (any(atoms(:)%is_qm == 1)) then
            filename = trim(base_filename) // '_qm_mm.txt'
            open(newunit=unit, file=filename, status='replace', action='write')
            write(unit, '(A)') '# QM/MM Region Marked Schlegel Diagram'
            write(unit, '(A)') '# QM atoms marked with *'
            write(unit, '(A)') '# ID  Element  QM?  X  Y  Z  Schlegel_X  Schlegel_Y'
            do i = 1, natoms
                call c_f_string(atoms(i)%element, element_str)
                if (atoms(i)%is_qm == 1) then
                    write(unit, '(I5,1X,A2,1X,A1,5F10.4)') &
                        atoms(i)%id, element_str, '*', &
                        atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                        atoms(i)%schlegel_x, atoms(i)%schlegel_y
                else
                    write(unit, '(I5,1X,A2,1X,A1,5F10.4)') &
                        atoms(i)%id, element_str, ' ', &
                        atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                        atoms(i)%schlegel_x, atoms(i)%schlegel_y
                end if
            end do
            close(unit)
        end if
        
    end subroutine export_schlegel_projections
    
    ! Main export function - C-callable
    subroutine export_all_formats(filename_c, natoms, nbonds, atoms, bonds, &
                                has_box, box_lo, box_hi, projection_point, view_angle) &
                                bind(C, name="export_all_formats")
        character(c_char), intent(in) :: filename_c(*)
        integer(c_int), intent(in), value :: natoms, nbonds, has_box
        type(export_atom_t), intent(in) :: atoms(*)  ! C assumed-size is OK here
        type(export_bond_t), intent(in) :: bonds(*)  ! C assumed-size is OK here
        real(c_double), intent(in) :: box_lo(3), box_hi(3)
        real(c_double), intent(in) :: projection_point(3)
        real(c_double), intent(in), value :: view_angle
        
        character(len=MAX_FILENAME_LEN) :: base_filename
        character(len=MAX_FILENAME_LEN) :: filename
        logical :: l_has_box
        type(export_atom_t), allocatable :: atoms_local(:)
        type(export_bond_t), allocatable :: bonds_local(:)
        integer :: i
        
        ! Convert C string to Fortran string
        call c_f_string(filename_c, base_filename)
        
        ! Create local arrays for modifications (now with explicit shape)
        allocate(atoms_local(natoms))
        allocate(bonds_local(nbonds))
        
        do i = 1, natoms
            atoms_local(i) = atoms(i)
        end do
        
        do i = 1, nbonds
            bonds_local(i) = bonds(i)
        end do
        
        l_has_box = (has_box /= 0)
        
        ! Calculate Schlegel projections
        call calculate_schlegel_projections(natoms, atoms_local, projection_point, view_angle)
        
        ! Export standard formats - now using local arrays (explicit shape)
        filename = trim(base_filename) // '.xyz'
        call export_xyz(filename, natoms, atoms_local, .true.)
        
        filename = trim(base_filename) // '.pdb'
        call export_pdb(filename, natoms, atoms_local, bonds_local, nbonds)
        
        filename = trim(base_filename) // '.sdf'
        call export_sdf(filename, natoms, nbonds, atoms_local, bonds_local)
        
        filename = trim(base_filename) // '.mol2'
        call export_mol2(filename, natoms, nbonds, atoms_local, bonds_local)
        
        filename = trim(base_filename) // '.cif'
        call export_cif(filename, natoms, atoms_local, l_has_box, box_lo, box_hi)
        
        filename = trim(base_filename) // '.gro'
        call export_gro(filename, natoms, atoms_local, l_has_box, box_lo, box_hi)
        
        ! Export Schlegel projections
        call export_schlegel_projections(base_filename, natoms, atoms_local, &
                                        projection_point, view_angle)
        
        deallocate(atoms_local, bonds_local)
        
    end subroutine export_all_formats
    
    ! Simplified export - just Schlegel projections
    subroutine export_schlegel_only(filename_c, natoms, atoms, projection_point, view_angle) &
                                   bind(C, name="export_schlegel_only")
        character(c_char), intent(in) :: filename_c(*)
        integer(c_int), intent(in) :: natoms
        type(export_atom_t), intent(in) :: atoms(*)
        real(c_double), intent(in) :: projection_point(3)
        real(c_double), intent(in) :: view_angle
        
        character(len=MAX_FILENAME_LEN) :: base_filename
        type(export_atom_t), allocatable :: atoms_local(:)
        integer :: i
        
        call c_f_string(filename_c, base_filename)
        
        allocate(atoms_local(natoms))
        do i = 1, natoms
            atoms_local(i) = atoms(i)
        end do
        
        call calculate_schlegel_projections(natoms, atoms_local, projection_point, view_angle)
        call export_schlegel_projections(base_filename, natoms, atoms_local, &
                                         projection_point, view_angle)
        
        deallocate(atoms_local)
        
    end subroutine export_schlegel_only
    
    ! Get version information
    subroutine get_version(version_str) bind(C, name="get_version")
        character(c_char), intent(out) :: version_str(*)
        call f_c_string('ADOBMD Molecular Exporter v2.0.0 with Schlegel projections', version_str)
    end subroutine get_version

end module molecular_export