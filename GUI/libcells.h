#ifndef LIBBCELLS_H
#define LIBBCELLS_H

#define DllExport __declspec( dllexport )

#ifdef __cplusplus
extern "C" {
#endif

void execute(int *, char *, int *,char *, int *);
void simulate_step(int *);
void terminate_run(int *);
void get_dimensions(int *,int *,int *, int *, double *, int *, bool *, double *);
void get_scene(int *, int *);
void get_summary(int *, int *, int *);
void get_concdata(int *, double *, double *);
void get_volprob(int *, double *, double *, double*);
void get_oxyprob(int *, double *, double*);
void CellML_load_file(char *, int, int *, int *);
void CellML_get_component_info(int, char *, int *);
void CellML_get_variable_name(int, int, char *);
void CellML_get_variable_value(int, int, double *);
void CellML_set_variable_value(int, int, double);

#ifdef __cplusplus
}
#endif

#endif // LIBBCELLS_H
