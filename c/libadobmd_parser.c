/* libadobmd_parser.c - ADOBMD parser in C
 * Compile: gcc -shared -fPIC -O3 -o libadobmd_parser.so libadobmd_parser.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define MAX_LINE_LEN 256
#define MAX_FILENAME_LEN 256

/* Atom structure (C-compatible) */
typedef struct {
    int id;
    int type_id;
    int molecule;
    char element[3];
    double x, y, z;
    double charge;
    int is_qm;
} atom_t;

/* Bond structure */
typedef struct {
    int id;
    int type_id;
    int atom1;
    int atom2;
    int order;
} bond_t;

/* Parser result */
typedef struct {
    int natoms;
    int nbonds;
    int nqm;
    double box_lo[3];
    double box_hi[3];
    int has_box;
    atom_t* atoms;
    bond_t* bonds;
    char filename[MAX_FILENAME_LEN];
} parse_result_t;

/* Helper: read line, skip comments and empty lines */
static int read_next_line(FILE* file, char* line, size_t size) {
    while (fgets(line, size, file)) {
        /* Skip empty lines */
        if (line[0] == '\n' || line[0] == '\r') continue;
        
        /* Skip comments */
        if (line[0] == '#') continue;
        
        return 1;
    }
    return 0;
}

/* Parse header */
static void parse_header(FILE* file, int* natoms, int* nbonds, int* natom_types,
                         bool* has_box, double* box_lo, double* box_hi) {
    char line[MAX_LINE_LEN];
    long pos;
    
    *natoms = 0;
    *nbonds = 0;
    *natom_types = 0;
    *has_box = false;
    box_lo[0] = box_lo[1] = box_lo[2] = 0.0;
    box_hi[0] = box_hi[1] = box_hi[2] = 0.0;
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '#') continue;
        
        if (strstr(line, "atoms")) {
            sscanf(line, "%d", natoms);
        }
        else if (strstr(line, "bonds")) {
            sscanf(line, "%d", nbonds);
        }
        else if (strstr(line, "atom types")) {
            sscanf(line, "%d", natom_types);
        }
        else if (strstr(line, "xlo xhi")) {
            sscanf(line, "%lf %lf", &box_lo[0], &box_hi[0]);
            *has_box = true;
        }
        else if (strstr(line, "ylo yhi")) {
            sscanf(line, "%lf %lf", &box_lo[1], &box_hi[1]);
        }
        else if (strstr(line, "zlo zhi")) {
            sscanf(line, "%lf %lf", &box_lo[2], &box_hi[2]);
        }
        else if (strstr(line, "Masses")) {
            break;
        }
    }
}

/* Parse type elements from Masses section */
static void parse_type_elements(FILE* file, int natom_types, char** type_elements) {
    char line[MAX_LINE_LEN];
    int type_id, i;
    double mass;
    char* comment;
    
    /* Find Masses section */
    rewind(file);
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "Masses")) break;
    }
    
    /* Skip blank lines */
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '#' && line[0] != '\n' && line[0] != '\r') break;
    }
    
    /* Read masses */
    for (i = 0; i < natom_types; i++) {
        if (i > 0) fgets(line, sizeof(line), file);
        
        sscanf(line, "%d %lf", &type_id, &mass);
        
        /* Extract element from comment */
        type_elements[type_id] = strdup("X");
        comment = strchr(line, '#');
        if (comment) {
            comment++; /* Skip # */
            while (*comment == ' ') comment++;
            if (*comment) {
                /* Take first 1-2 letters as element */
                if (comment[1] >= 'a' && comment[1] <= 'z') {
                    type_elements[type_id] = malloc(3);
                    type_elements[type_id][0] = comment[0];
                    type_elements[type_id][1] = comment[1];
                    type_elements[type_id][2] = '\0';
                } else {
                    type_elements[type_id] = malloc(2);
                    type_elements[type_id][0] = comment[0];
                    type_elements[type_id][1] = '\0';
                }
            }
        }
    }
}

/* Parse atoms section */
static void parse_atoms(FILE* file, int natoms, atom_t* atoms,
                        char** type_elements, int* qm_indices, int* nqm) {
    char line[MAX_LINE_LEN];
    int i, atom_id, molecule, type_id;
    double x, y, z, charge;
    
    /* Find Atoms section */
    rewind(file);
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "Atoms")) break;
    }
    
    /* Skip blank lines */
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '#' && line[0] != '\n' && line[0] != '\r') break;
    }
    
    *nqm = 0;
    for (i = 0; i < natoms; i++) {
        if (i > 0) fgets(line, sizeof(line), file);
        
        sscanf(line, "%d %d %d %lf %lf %lf %lf",
               &atom_id, &molecule, &type_id, &charge, &x, &y, &z);
        
        atoms[i].id = atom_id;
        atoms[i].type_id = type_id;
        atoms[i].molecule = molecule;
        atoms[i].x = x;
        atoms[i].y = y;
        atoms[i].z = z;
        atoms[i].charge = charge;
        
        /* Set element from type */
        if (type_elements[type_id]) {
            strncpy(atoms[i].element, type_elements[type_id], 2);
            atoms[i].element[2] = '\0';
        } else {
            strcpy(atoms[i].element, "X");
        }
        
        /* Check for QM region */
        atoms[i].is_qm = 0;
        if (strstr(line, "QM") || strstr(line, " Q ")) {
            atoms[i].is_qm = 1;
            qm_indices[(*nqm)++] = atom_id;
        }
    }
}

/* Parse bonds section */
static void parse_bonds(FILE* file, int nbonds, bond_t* bonds) {
    char line[MAX_LINE_LEN];
    int i;
    
    if (nbonds == 0) return;
    
    /* Find Bonds section */
    rewind(file);
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "Bonds")) break;
    }
    
    /* Skip blank lines */
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '#' && line[0] != '\n' && line[0] != '\r') break;
    }
    
    for (i = 0; i < nbonds; i++) {
        if (i > 0) fgets(line, sizeof(line), file);
        
        sscanf(line, "%d %d %d %d",
               &bonds[i].id, &bonds[i].type_id,
               &bonds[i].atom1, &bonds[i].atom2);
        
        bonds[i].order = bonds[i].type_id;
    }
}

/* Main parsing function - C-callable */
void parse_file(const char* filename_c, parse_result_t* result_c) {
    FILE* file;
    int natoms = 0, nbonds = 0, natom_types = 0, nqm = 0;
    double box_lo[3], box_hi[3];
    bool has_box;
    atom_t* atoms = NULL;
    bond_t* bonds = NULL;
    char** type_elements = NULL;
    int* qm_indices = NULL;
    int i;
    
    /* Initialize result */
    result_c->natoms = 0;
    result_c->nbonds = 0;
    result_c->nqm = 0;
    result_c->atoms = NULL;
    result_c->bonds = NULL;
    strncpy(result_c->filename, filename_c, MAX_FILENAME_LEN - 1);
    result_c->filename[MAX_FILENAME_LEN - 1] = '\0';
    
    /* Open file */
    file = fopen(filename_c, "r");
    if (!file) {
        result_c->natoms = -1;
        return;
    }
    
    /* Parse header */
    parse_header(file, &natoms, &nbonds, &natom_types, &has_box, box_lo, box_hi);
    
    /* Allocate memory */
    atoms = (atom_t*)calloc(natoms, sizeof(atom_t));
    bonds = (bond_t*)calloc(nbonds, sizeof(bond_t));
    type_elements = (char**)calloc(natom_types + 1, sizeof(char*));
    qm_indices = (int*)calloc(natoms, sizeof(int));
    
    /* Parse type elements */
    parse_type_elements(file, natom_types, type_elements);
    
    /* Parse atoms */
    parse_atoms(file, natoms, atoms, type_elements, qm_indices, &nqm);
    
    /* Parse bonds */
    if (nbonds > 0) {
        parse_bonds(file, nbonds, bonds);
    }
    
    fclose(file);
    
    /* Fill result structure */
    result_c->natoms = natoms;
    result_c->nbonds = nbonds;
    result_c->nqm = nqm;
    result_c->box_lo[0] = box_lo[0];
    result_c->box_lo[1] = box_lo[1];
    result_c->box_lo[2] = box_lo[2];
    result_c->box_hi[0] = box_hi[0];
    result_c->box_hi[1] = box_hi[1];
    result_c->box_hi[2] = box_hi[2];
    result_c->has_box = has_box ? 1 : 0;
    result_c->atoms = atoms;
    result_c->bonds = bonds;
    
    /* Free temporary arrays */
    for (i = 0; i <= natom_types; i++) {
        if (type_elements[i]) free(type_elements[i]);
    }
    free(type_elements);
    free(qm_indices);
}

/* Free memory function */
void free_result(parse_result_t* result_c) {
    if (result_c->atoms) {
        free(result_c->atoms);
        result_c->atoms = NULL;
    }
    
    if (result_c->bonds) {
        free(result_c->bonds);
        result_c->bonds = NULL;
    }
    
    result_c->natoms = 0;
    result_c->nbonds = 0;
    result_c->nqm = 0;
}

/* Get version information */
void get_version(char* version_str) {
    strcpy(version_str, "ADOBMD C Parser v1.0.0");
}