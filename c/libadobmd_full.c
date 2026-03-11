/* libadobmd_full.c - Complete ADOBMD parser + exporter with Schlegel projections
 * Compile: gcc -shared -fPIC -O3 -lm -o libadobmd_full.so libadobmd_full.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define MAX_LINE_LEN 256
#define MAX_FILENAME_LEN 256
#define MAX_ELEMENT_LEN 3
#define PI 3.14159265358979323846
#define DEG_TO_RAD (PI / 180.0)

/* Atom structure with Schlegel projections */
typedef struct {
    int id;
    int type_id;
    int molecule;
    char element[4];  /* +1 for null terminator */
    double x, y, z;
    double charge;
    int is_qm;
    /* Projections for Schlegel diagrams */
    double xy_x, xy_y;      /* XY plane (top view) */
    double xz_x, xz_z;      /* XZ plane (side view) */
    double yz_y, yz_z;      /* YZ plane (side view) */
    double schlegel_x, schlegel_y;  /* Main Schlegel view */
} export_atom_t;

/* Bond structure */
typedef struct {
    int id;
    int type_id;
    int atom1;
    int atom2;
    int order;
} export_bond_t;

/* Main data structure */
typedef struct {
    int natoms;
    int nbonds;
    int nqm;
    double box_lo[3];
    double box_hi[3];
    int has_box;
    export_atom_t* atoms;
    export_bond_t* bonds;
    char filename[MAX_FILENAME_LEN];
    /* Schlegel parameters */
    double projection_point[3];
    double view_angle;
    double scale_factor;
} molecular_data_t;

/* ============ PARSING FUNCTIONS ============ */

static int read_next_line(FILE* file, char* line) {
    while (fgets(line, MAX_LINE_LEN, file)) {
        char* p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r') continue;
        return 1;
    }
    return 0;
}

static void parse_header(FILE* file, int* natoms, int* nbonds, int* natom_types,
                         bool* has_box, double* box_lo, double* box_hi) {
    char line[MAX_LINE_LEN];
    
    *natoms = 0; *nbonds = 0; *natom_types = 0;
    *has_box = false;
    box_lo[0] = box_lo[1] = box_lo[2] = 0.0;
    box_hi[0] = box_hi[1] = box_hi[2] = 0.0;
    
    while (fgets(line, MAX_LINE_LEN, file)) {
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

static void parse_type_elements(FILE* file, int natom_types, char** type_elements) {
    char line[MAX_LINE_LEN];
    int type_id, i;
    double mass;
    char* comment;
    
    rewind(file);
    while (fgets(line, MAX_LINE_LEN, file)) {
        if (strstr(line, "Masses")) break;
    }
    
    for (i = 0; i < natom_types; i++) {
        if (!fgets(line, MAX_LINE_LEN, file)) break;
        
        sscanf(line, "%d %lf", &type_id, &mass);
        
        /* Default to X */
        type_elements[type_id] = strdup("X");
        
        /* Extract element from comment */
        comment = strchr(line, '#');
        if (comment) {
            comment++;
            while (*comment == ' ') comment++;
            if (*comment) {
                if (comment[1] >= 'a' && comment[1] <= 'z') {
                    /* Two-letter element */
                    type_elements[type_id] = malloc(3);
                    type_elements[type_id][0] = comment[0];
                    type_elements[type_id][1] = comment[1];
                    type_elements[type_id][2] = '\0';
                } else {
                    /* One-letter element */
                    type_elements[type_id] = malloc(2);
                    type_elements[type_id][0] = comment[0];
                    type_elements[type_id][1] = '\0';
                }
            }
        }
    }
}

static void parse_atoms(FILE* file, int natoms, export_atom_t* atoms,
                        char** type_elements, int** qm_indices, int* nqm) {
    char line[MAX_LINE_LEN];
    int i;
    
    rewind(file);
    while (fgets(line, MAX_LINE_LEN, file)) {
        if (strstr(line, "Atoms")) break;
    }
    
    *nqm = 0;
    for (i = 0; i < natoms; i++) {
        if (!fgets(line, MAX_LINE_LEN, file)) break;
        
        export_atom_t* a = &atoms[i];
        sscanf(line, "%d %d %d %lf %lf %lf %lf",
               &a->id, &a->molecule, &a->type_id, 
               &a->charge, &a->x, &a->y, &a->z);
        
        /* Set element */
        if (type_elements[a->type_id]) {
            strcpy(a->element, type_elements[a->type_id]);
        } else {
            strcpy(a->element, "X");
        }
        
        /* Check QM */
        a->is_qm = (strstr(line, "QM") || strstr(line, " Q ")) ? 1 : 0;
        if (a->is_qm) {
            (*qm_indices)[*nqm] = a->id;
            (*nqm)++;
        }
    }
}

static void parse_bonds(FILE* file, int nbonds, export_bond_t* bonds) {
    char line[MAX_LINE_LEN];
    int i;
    
    if (nbonds == 0) return;
    
    rewind(file);
    while (fgets(line, MAX_LINE_LEN, file)) {
        if (strstr(line, "Bonds")) break;
    }
    
    for (i = 0; i < nbonds; i++) {
        if (!fgets(line, MAX_LINE_LEN, file)) break;
        
        export_bond_t* b = &bonds[i];
        sscanf(line, "%d %d %d %d",
               &b->id, &b->type_id, &b->atom1, &b->atom2);
        b->order = b->type_id;
    }
}

/* ============ SCHLEGEL PROJECTIONS ============ */

static void perspective_projection(double point[3], double proj_point[3], 
                                   double scale, double proj_out[2]) {
    double z_dist, factor;
    
    z_dist = proj_point[2] - point[2];
    if (fabs(z_dist) < 1e-10) {
        proj_out[0] = point[0] * scale;
        proj_out[1] = point[1] * scale;
    } else {
        factor = (proj_point[2] - point[2]) / z_dist * scale;
        proj_out[0] = point[0] * factor;
        proj_out[1] = point[1] * factor;
    }
}

static void calculate_schlegel_projections(molecular_data_t* data) {
    int i;
    double point[3], proj_out[2];
    
    for (i = 0; i < data->natoms; i++) {
        export_atom_t* a = &data->atoms[i];
        point[0] = a->x; point[1] = a->y; point[2] = a->z;
        
        /* Basic projections */
        a->xy_x = point[0];
        a->xy_y = point[1];
        
        a->xz_x = point[0];
        a->xz_z = point[2];
        
        a->yz_y = point[1];
        a->yz_z = point[2];
        
        /* Perspective projection for main Schlegel view */
        perspective_projection(point, data->projection_point, 1.0, proj_out);
        a->schlegel_x = proj_out[0];
        a->schlegel_y = proj_out[1];
    }
}

/* ============ EXPORT FUNCTIONS ============ */

static void export_xyz(const char* base_filename, molecular_data_t* data) {
    char filename[MAX_FILENAME_LEN];
    FILE* f;
    int i;
    
    snprintf(filename, MAX_FILENAME_LEN, "%s.xyz", base_filename);
    f = fopen(filename, "w");
    if (!f) return;
    
    fprintf(f, "%d\n", data->natoms);
    fprintf(f, "ADOBMD Export with QM/MM regions\n");
    
    for (i = 0; i < data->natoms; i++) {
        export_atom_t* a = &data->atoms[i];
        fprintf(f, "%-2s %15.8f %15.8f %15.8f\n", 
                a->element, a->x, a->y, a->z);
    }
    
    fclose(f);
}

static void export_schlegel_projections(const char* base_filename, molecular_data_t* data) {
    char filename[MAX_FILENAME_LEN];
    FILE* f;
    int i;
    
    /* XY plane */
    snprintf(filename, MAX_FILENAME_LEN, "%s_xy.txt", base_filename);
    f = fopen(filename, "w");
    fprintf(f, "# XY Plane Projection (top view)\n");
    fprintf(f, "# Projection point: %f %f %f\n", 
            data->projection_point[0], data->projection_point[1], data->projection_point[2]);
    fprintf(f, "# ID Element X Y Z Proj_X Proj_Y\n");
    for (i = 0; i < data->natoms; i++) {
        export_atom_t* a = &data->atoms[i];
        fprintf(f, "%5d %-2s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                a->id, a->element, a->x, a->y, a->z, a->xy_x, a->xy_y);
    }
    fclose(f);
    
    /* XZ plane */
    snprintf(filename, MAX_FILENAME_LEN, "%s_xz.txt", base_filename);
    f = fopen(filename, "w");
    fprintf(f, "# XZ Plane Projection (side view)\n");
    fprintf(f, "# ID Element X Y Z Proj_X Proj_Z\n");
    for (i = 0; i < data->natoms; i++) {
        export_atom_t* a = &data->atoms[i];
        fprintf(f, "%5d %-2s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                a->id, a->element, a->x, a->y, a->z, a->xz_x, a->xz_z);
    }
    fclose(f);
    
    /* YZ plane */
    snprintf(filename, MAX_FILENAME_LEN, "%s_yz.txt", base_filename);
    f = fopen(filename, "w");
    fprintf(f, "# YZ Plane Projection (side view)\n");
    fprintf(f, "# ID Element X Y Z Proj_Y Proj_Z\n");
    for (i = 0; i < data->natoms; i++) {
        export_atom_t* a = &data->atoms[i];
        fprintf(f, "%5d %-2s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                a->id, a->element, a->x, a->y, a->z, a->yz_y, a->yz_z);
    }
    fclose(f);
    
    /* Main Schlegel view */
    snprintf(filename, MAX_FILENAME_LEN, "%s_schlegel.txt", base_filename);
    f = fopen(filename, "w");
    fprintf(f, "# Schlegel Diagram (perspective projection)\n");
    fprintf(f, "# ID Element X Y Z Schlegel_X Schlegel_Y\n");
    for (i = 0; i < data->natoms; i++) {
        export_atom_t* a = &data->atoms[i];
        fprintf(f, "%5d %-2s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                a->id, a->element, a->x, a->y, a->z, a->schlegel_x, a->schlegel_y);
    }
    fclose(f);
    
    /* QM/MM colored version */
    if (data->nqm > 0) {
        snprintf(filename, MAX_FILENAME_LEN, "%s_qm_mm.txt", base_filename);
        f = fopen(filename, "w");
        fprintf(f, "# QM/MM Region Marked Schlegel Diagram\n");
        fprintf(f, "# QM atoms marked with *\n");
        fprintf(f, "# ID Element QM? X Y Z Schlegel_X Schlegel_Y\n");
        for (i = 0; i < data->natoms; i++) {
            export_atom_t* a = &data->atoms[i];
            fprintf(f, "%5d %-2s %c %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                    a->id, a->element, a->is_qm ? '*' : ' ',
                    a->x, a->y, a->z, a->schlegel_x, a->schlegel_y);
        }
        fclose(f);
    }
}

/* ============ C INTERFACE ============ */

/* Parse file and return data */
molecular_data_t* parse_adobmd_file(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) return NULL;
    
    /* Allocate data structure */
    molecular_data_t* data = calloc(1, sizeof(molecular_data_t));
    strncpy(data->filename, filename, MAX_FILENAME_LEN - 1);
    
    /* Default Schlegel parameters */
    data->projection_point[0] = 0.0;
    data->projection_point[1] = 0.0;
    data->projection_point[2] = 10.0;  /* Project from above */
    data->view_angle = 45.0;
    data->scale_factor = 1.0;
    
    /* Parse header */
    int natom_types;
    bool has_box;
    parse_header(file, &data->natoms, &data->nbonds, &natom_types, 
                 &has_box, data->box_lo, data->box_hi);
    data->has_box = has_box ? 1 : 0;
    
    if (data->natoms <= 0) {
        fclose(file);
        free(data);
        return NULL;
    }
    
    /* Allocate arrays */
    data->atoms = calloc(data->natoms, sizeof(export_atom_t));
    data->bonds = calloc(data->nbonds, sizeof(export_bond_t));
    char** type_elements = calloc(natom_types + 1, sizeof(char*));
    int* qm_indices = calloc(data->natoms, sizeof(int));
    
    /* Parse sections */
    parse_type_elements(file, natom_types, type_elements);
    parse_atoms(file, data->natoms, data->atoms, type_elements, &qm_indices, &data->nqm);
    parse_bonds(file, data->nbonds, data->bonds);
    
    fclose(file);
    
    /* Free temporary arrays */
    for (int i = 0; i <= natom_types; i++) {
        if (type_elements[i]) free(type_elements[i]);
    }
    free(type_elements);
    free(qm_indices);
    
    /* Calculate Schlegel projections */
    calculate_schlegel_projections(data);
    
    return data;
}

/* Export all formats */
void export_all_formats(molecular_data_t* data, const char* base_filename) {
    if (!data || !base_filename) return;
    
    export_xyz(base_filename, data);
    export_schlegel_projections(base_filename, data);
    
    /* Add more export formats here (PDB, SDF, MOL2, CIF, GRO) */
}

/* Set Schlegel projection parameters */
void set_schlegel_params(molecular_data_t* data, 
                         double proj_x, double proj_y, double proj_z,
                         double view_angle, double scale) {
    if (!data) return;
    
    data->projection_point[0] = proj_x;
    data->projection_point[1] = proj_y;
    data->projection_point[2] = proj_z;
    data->view_angle = view_angle;
    data->scale_factor = scale;
    
    /* Recalculate projections */
    calculate_schlegel_projections(data);
}

/* Free data */
void free_molecular_data(molecular_data_t* data) {
    if (!data) return;
    
    if (data->atoms) free(data->atoms);
    if (data->bonds) free(data->bonds);
    free(data);
}

/* Get atom data for Python */
void get_atom_data(molecular_data_t* data, int index, 
                   int* id, char* element, double* x, double* y, double* z,
                   double* charge, int* is_qm,
                   double* schlegel_x, double* schlegel_y) {
    if (!data || index < 0 || index >= data->natoms) return;
    
    export_atom_t* a = &data->atoms[index];
    *id = a->id;
    strcpy(element, a->element);
    *x = a->x; *y = a->y; *z = a->z;
    *charge = a->charge;
    *is_qm = a->is_qm;
    *schlegel_x = a->schlegel_x;
    *schlegel_y = a->schlegel_y;
}

/* Get version */
const char* get_version(void) {
    return "ADOBMD Full Parser/Exporter v1.0.0 with Schlegel projections";
}