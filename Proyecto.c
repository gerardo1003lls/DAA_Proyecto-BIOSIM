#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sqlite3.h"

#define NUM_TERRITORIOS 20
#define SIN_CONEXION 0.0
#define MAX_INDIVIDUOS 150

typedef struct Contacto{
    int u_individuo;
    int v_individuo;
    float prob_contagio;
    struct Contacto *sgt;
} Contacto;

typedef struct Individuo{
    int ID;
    char Nombre[50];
    int Territorio_ID;
    float Riesgo_inicial;
    int Grado_inicial;
    
    int Infectado;
    int t_infeccion;
    int Cepa_ID;
    int Recuperado;
    int Fallecido;
    
    struct Contacto *contactos;
}Individuo;

typedef struct Territorio{
    int ID;
    char Nombre[50];
    int M;
    
    struct Individuo *individuos[MAX_INDIVIDUOS];
    int num_individuos;
}Territorio;

typedef struct ConexionTerritorio{
    int u_territorio;
    int v_territorio;
    float peso_proximidad;
}ConexionTerritorio;

typedef struct Mapa{
    Territorio territorios[NUM_TERRITORIOS];
    float matrix[NUM_TERRITORIOS][NUM_TERRITORIOS];
    int num_territorios;
    int num_conexiones;
}Mapa;

enum TerritoriosIdx {
    CHINA = 0, JAPON, KOREA, TURQUIA, ARABIA, RUSIA, GRECIA, CROACIA, HUNGRIA, POLONIA,
    FINLANDIA, SUECIA, DINAMARCA, ALEMANIA, FRANCIA, ESPANA, PORTUGAL, ITALIA, EUA, REINO_UNIDO
};

int IDs = 0;

void MENU();
void CrearTerritorio(Territorio *t, int id, const char *nom, int cap);
void CrearConexiones(Mapa *grafo);
float Azar(float lim_inf, float lim_sup);
void AgregarConexion(Mapa *grafo, int t1, int t2, float peso);
void InicializarGrafo(Mapa *grafo);
void AgregarIndividuo(Territorio *territorio, Individuo *individuo);
void CrearIndividuos(Mapa *grafo, int territorio_id, sqlite3 *db);
void ConsultaSQL(int territorio_id, char *buffer);

int main(int argc, char const *argv[]){
    srand(42);
    
    sqlite3 *db;
    int rc = sqlite3_open("Nombres.db", &db);
    
    Mapa mundo;
    InicializarGrafo(&mundo);
    CrearConexiones(&mundo);
    
    // Creamos los individuos de los territorios
    for(int i = 0; i < NUM_TERRITORIOS; i++){
        CrearIndividuos(&mundo, i, db);
    }
    
    MENU();
    
    sqlite3_close(db);
    getchar();
    return 0;
}

void CrearTerritorio(Territorio *t, int id, const char *nom, int cap){
    t->ID = id;
    strcpy(t->Nombre, nom);
    t->M = cap;
    t->num_individuos = 0;
    
    for(int i = 0; i < MAX_INDIVIDUOS; i++) {
        t->individuos[i] = NULL;
    }
}

void CrearConexiones(Mapa *grafo){
    AgregarConexion(grafo, CHINA, KOREA, 0.85);
    AgregarConexion(grafo, CHINA, JAPON, 0.75);
    AgregarConexion(grafo, CHINA, RUSIA, 0.80);
    AgregarConexion(grafo, KOREA, JAPON, 0.70);
    AgregarConexion(grafo, JAPON, RUSIA, 0.65);
    AgregarConexion(grafo, RUSIA, TURQUIA, 0.75);
    
    AgregarConexion(grafo, TURQUIA, ARABIA, 0.80);
    AgregarConexion(grafo, TURQUIA, GRECIA, 0.90);
    AgregarConexion(grafo, ARABIA, GRECIA, 0.60);
    
    AgregarConexion(grafo, RUSIA, POLONIA, 0.85);
    AgregarConexion(grafo, RUSIA, FINLANDIA, 0.88);
    AgregarConexion(grafo, POLONIA, HUNGRIA, 0.82);
    AgregarConexion(grafo, POLONIA, ALEMANIA, 0.90);
    AgregarConexion(grafo, HUNGRIA, CROACIA, 0.85);
    AgregarConexion(grafo, HUNGRIA, ITALIA, 0.75);
    
    AgregarConexion(grafo, FINLANDIA, SUECIA, 0.92);
    AgregarConexion(grafo, SUECIA, DINAMARCA, 0.88);
    AgregarConexion(grafo, DINAMARCA, ALEMANIA, 0.90);
    AgregarConexion(grafo, SUECIA, POLONIA, 0.70);
    
    AgregarConexion(grafo, ALEMANIA, FRANCIA, 0.92);
    AgregarConexion(grafo, ALEMANIA, ITALIA, 0.80);
    AgregarConexion(grafo, FRANCIA, ESPANA, 0.88);
    AgregarConexion(grafo, FRANCIA, ITALIA, 0.85);
    AgregarConexion(grafo, ESPANA, PORTUGAL, 0.95);
    AgregarConexion(grafo, ITALIA, GRECIA, 0.72);
    AgregarConexion(grafo, ITALIA, CROACIA, 0.83);
    
    AgregarConexion(grafo, REINO_UNIDO, FRANCIA, 0.87);
    AgregarConexion(grafo, REINO_UNIDO, ESPANA, 0.65);
    
    AgregarConexion(grafo, EUA, REINO_UNIDO, 0.78);
    AgregarConexion(grafo, EUA, JAPON, 0.70);
}

void AgregarConexion(Mapa *grafo, int t1, int t2, float peso){
    grafo->matrix[t1][t2] = peso;
    grafo->matrix[t2][t1] = peso;
    grafo->num_conexiones++;
}

void InicializarGrafo(Mapa *grafo){
    grafo->num_territorios = NUM_TERRITORIOS;
    grafo->num_conexiones = 0;
    
    CrearTerritorio(&grafo->territorios[CHINA], 0, "China", 150);
    CrearTerritorio(&grafo->territorios[JAPON], 1, "Japon", 120);
    CrearTerritorio(&grafo->territorios[KOREA], 2, "Korea", 100);
    CrearTerritorio(&grafo->territorios[TURQUIA], 3, "Turquia", 110);
    CrearTerritorio(&grafo->territorios[ARABIA], 4, "Arabia", 90);
    CrearTerritorio(&grafo->territorios[RUSIA], 5, "Rusia", 180);
    CrearTerritorio(&grafo->territorios[GRECIA], 6, "Grecia", 80);
    CrearTerritorio(&grafo->territorios[CROACIA], 7, "Croacia", 70);
    CrearTerritorio(&grafo->territorios[HUNGRIA], 8, "Hungria", 85);
    CrearTerritorio(&grafo->territorios[POLONIA], 9, "Polonia", 95);
    CrearTerritorio(&grafo->territorios[FINLANDIA], 10, "Finlandia", 75);
    CrearTerritorio(&grafo->territorios[SUECIA], 11, "Suecia", 90);
    CrearTerritorio(&grafo->territorios[DINAMARCA], 12, "Dinamarca", 80);
    CrearTerritorio(&grafo->territorios[ALEMANIA], 13, "Alemania", 140);
    CrearTerritorio(&grafo->territorios[FRANCIA], 14, "Francia", 130);
    CrearTerritorio(&grafo->territorios[ESPANA], 15, "Espana", 115);
    CrearTerritorio(&grafo->territorios[PORTUGAL], 16, "Portugal", 85);
    CrearTerritorio(&grafo->territorios[ITALIA], 17, "Italia", 120);
    CrearTerritorio(&grafo->territorios[EUA], 18, "EUA", 200);
    CrearTerritorio(&grafo->territorios[REINO_UNIDO], 19, "Reino Unido", 125);
    
    for(int i = 0; i < NUM_TERRITORIOS; i++) {
        for(int j = 0; j < NUM_TERRITORIOS; j++) {
            grafo->matrix[i][j] = SIN_CONEXION;
        }
    }
}

void MENU(){
    printf("\n============ BioSim ============");
    printf("\n[1] Inicializacion y analisis de datos");
    printf("\n[2] Deteccion de brotes");
    printf("\n[3] Propagacion temporal");
    printf("\n[4] Minimizacion del riesgo total ");
    printf("\n[5] Identificacion de rutas criticas");
    printf("\n[6] Calculo de rutas optimas de contencion");
    printf("\n[7] Clustering de cepas similares");
    printf("\n[8] Almacenamiento eficiente y consulta rapida");
    printf("\n[0] Salir");
    printf("\n================================\n");
}

void AgregarIndividuo(Territorio *territorio, Individuo *individuo){
    territorio->individuos[territorio->num_individuos] = individuo;
    territorio->num_individuos++;
    individuo->Territorio_ID = territorio->ID;
}

void CrearIndividuos(Mapa *grafo, int territorio_id, sqlite3 *db){
    sqlite3_stmt *stmt;
    char sql[200];
    ConsultaSQL(territorio_id, sql);
    
    int rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    
    Territorio *territorio = &grafo->territorios[territorio_id];
    
    for(int i = 0; i < territorio->M && sqlite3_step(stmt) == SQLITE_ROW; i++){
        Individuo *P = (Individuo*)malloc(sizeof(Individuo));
        
        const unsigned char *nombre = sqlite3_column_text(stmt, 0);
        
        P->ID = IDs;
        strcpy(P->Nombre, (const char*)nombre);
        P->Territorio_ID = territorio_id;
        P->Riesgo_inicial = Azar(0.0, 1.0);
        P->Grado_inicial = (int)Azar(1, 10);
        P->Infectado = 0;
        P->t_infeccion = -1;
        P->Recuperado = 0;
        P->Fallecido = 0;
        P->Cepa_ID = -1;
        P->contactos = NULL;

        AgregarIndividuo(territorio, P);
        IDs++;
    }
    
    sqlite3_finalize(stmt);
}

void ConsultaSQL(int territorio_id, char *buffer){
    strcpy(buffer, "SELECT name FROM firstnames WHERE ");
    
    char territorio[30];
    switch(territorio_id){
        case 0: strcpy(territorio, "China"); break;
        case 1: strcpy(territorio, "Japan"); break;
        case 2: strcpy(territorio, "Korea"); break;
        case 3: strcpy(territorio, "Turkey"); break;
        case 4: strcpy(territorio, "[Arabia/Persia]"); break;
        case 5: strcpy(territorio, "Russia"); break;
        case 6: strcpy(territorio, "Greece"); break;
        case 7: strcpy(territorio, "Croatia"); break;
        case 8: strcpy(territorio, "Hungary"); break;
        case 9: strcpy(territorio, "Poland"); break;
        case 10: strcpy(territorio, "Finland"); break;
        case 11: strcpy(territorio, "Sweden"); break;
        case 12: strcpy(territorio, "Denmark"); break;
        case 13: strcpy(territorio, "Germany"); break;
        case 14: strcpy(territorio, "France"); break;
        case 15: strcpy(territorio, "Spain"); break;
        case 16: strcpy(territorio, "Portugal"); break;
        case 17: strcpy(territorio, "Italy"); break;
        case 18: strcpy(territorio, "[U.S.A.]"); break;
        case 19: strcpy(territorio, "[Great Britain]"); break;
        default: strcpy(territorio, "China"); break;
    }
    
    strcat(buffer, territorio);
    strcat(buffer, " IS NOT NULL LIMIT ");
    
    char limit[10];
    sprintf(limit, "%d", MAX_INDIVIDUOS);
    strcat(buffer, limit);
    strcat(buffer, ";");
}

float Azar(float lim_inf, float lim_sup){ 
    float num = rand() / (float) RAND_MAX;
    num = lim_inf + num * (lim_sup - lim_inf);
    return num;
}