//biosim - simulador de propagacion y contencion de epidemias
//proyecto final de algoritmos
//compilar: gcc Proyecto.c sqlite3.c -o BioSim -lm

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sqlite3.h"

//constantes del sistema
#define NUM_TERRITORIOS 20      //numero de paises/territorios
#define SIN_CONEXION 0.0        //valor para indicar que no hay conexion
#define MAX_INDIVIDUOS 150      //maximo de personas por territorio
#define NUM_CEPAS 50            //numero de variantes del virus
#define HASH_SIZE 2053          //tamanio de la tabla hash (numero primo)
#define ALPHABET_SIZE 26        //letras del alfabeto para el trie

//estructura para representar un contacto entre dos individuos
typedef struct Contacto{
    int u_individuo;            //id del primer individuo
    int v_individuo;            //id del segundo individuo
    float prob_contagio;        //probabilidad de contagio entre ellos
    struct Contacto *sgt;       //siguiente contacto en la lista
} Contacto;

//estructura para representar una persona en la simulacion
typedef struct Individuo{
    int ID;                     //identificador unico
    char Nombre[50];            //nombre de la persona
    int Territorio_ID;          //territorio donde vive
    float Riesgo_inicial;       //nivel de riesgo base
    int Grado_inicial;          //numero de contactos
    
    int Infectado;              //1 si esta infectado, 0 si no
    int t_infeccion;            //dia en que se infecto
    int Cepa_ID;                //cepa con la que esta infectado
    int Recuperado;             //1 si ya se recupero
    int Fallecido;              //1 si fallecio
    
    struct Contacto *contactos; //lista de contactos
}Individuo;

//estructura para representar una cepa del virus
typedef struct Cepa{
    int ID;                     //identificador de la cepa
    char Nombre[30];            //nombre de la variante
    float Tasa_contagio;        //probabilidad de contagiar (beta)
    float Tasa_mortalidad;      //probabilidad de muerte
    int Tiempo_incubacion;      //dias hasta mostrar sintomas
    int Tiempo_recuperacion;    //dias hasta recuperarse (gamma)
    char Sintomas[100];         //descripcion de sintomas
} Cepa;

//nodo para el heap minimo usado en dijkstra y prim
typedef struct NodoHeap{
    int vertice;                //indice del territorio
    float prioridad;            //valor de prioridad (menor = mas urgente)
} NodoHeap;

//estructura del heap minimo para algoritmos de grafos
typedef struct MinHeap{
    NodoHeap *nodos;            //arreglo de nodos
    int *posiciones;            //posicion de cada vertice en el heap
    int capacidad;              //tamanio maximo
    int size;                   //elementos actuales
} MinHeap;

//nodo para la tabla hash de individuos
typedef struct NodoHash{
    int ID;                         //id del individuo
    Individuo *individuo;           //puntero al individuo
    struct NodoHash *siguiente;     //siguiente nodo (para colisiones)
} NodoHash;

//tabla hash para busqueda de individuos en o(1)
typedef struct HashTable{
    NodoHash *tabla[HASH_SIZE];     //arreglo de punteros a nodos
    int num_elementos;              //cantidad de elementos almacenados
} HashTable;

//nodo del trie para clustering de cepas
typedef struct NodoTrie{
    struct NodoTrie *hijos[ALPHABET_SIZE];  //26 hijos (a-z)
    int es_final;                   //1 si termina una palabra aqui
    int cepa_id;                    //id de la cepa si es final
} NodoTrie;

//estructura principal del trie
typedef struct Trie{
    NodoTrie *raiz;                 //nodo raiz
    int num_palabras;               //cantidad de palabras insertadas
} Trie;

//estructura para los 10 pacientes infectados iniciales
typedef struct Semilla{
    int individuo_id;               //id del paciente cero
    int t0;                         //dia de inicio de infeccion
    int cepa_id;                    //cepa con la que se infecto
} Semilla;

//nodo para la tabla hash de cepas
typedef struct NodoHashCepa{
    int ID;                         //id de la cepa
    Cepa *cepa;                     //puntero a la cepa
    struct NodoHashCepa *siguiente; //siguiente nodo
} NodoHashCepa;

//tabla hash para busqueda de cepas en o(1)
typedef struct HashTableCepas{
    NodoHashCepa *tabla[NUM_CEPAS * 2];  //tamanio doble para evitar colisiones
    int num_elementos;
} HashTableCepas;

//estructura union-find para detectar componentes conexas
typedef struct UnionFind{
    int *padre;                     //arreglo de padres
    int *rango;                     //arreglo de rangos
    int n;                          //numero de elementos
} UnionFind;

//estructura para representar un territorio/pais
typedef struct Territorio{
    int ID;                         //identificador del territorio
    char Nombre[50];                //nombre del pais
    int M;                          //capacidad maxima de individuos
    
    struct Individuo *individuos[MAX_INDIVIDUOS];  //arreglo de habitantes
    int num_individuos;             //cantidad actual de habitantes
}Territorio;

//estructura para conexiones entre territorios
typedef struct ConexionTerritorio{
    int u_territorio;               //primer territorio
    int v_territorio;               //segundo territorio
    float peso_proximidad;          //peso de la conexion
}ConexionTerritorio;

//estructura principal que contiene todo el sistema
typedef struct Mapa{
    Territorio territorios[NUM_TERRITORIOS];  //arreglo de territorios
    float matrix[NUM_TERRITORIOS][NUM_TERRITORIOS];  //matriz de adyacencia
    int num_territorios;            //cantidad de territorios
    int num_conexiones;             //cantidad de conexiones
    
    Cepa cepas[NUM_CEPAS];          //arreglo de cepas
    int num_cepas;                  //cantidad de cepas
    
    HashTable *hash_individuos;     //hash para buscar individuos
    HashTableCepas *hash_cepas;     //hash para buscar cepas
    Trie *trie_cepas;               //trie para clustering de cepas
    
    Semilla semillas[10];           //10 semillas iniciales
    int num_semillas;               //cantidad de semillas
} Mapa;

//estructura de cola para bfs
typedef struct Cola{
    int *items;                     //arreglo de elementos
    int frente;                     //indice del frente
    int final;                      //indice del final
    int capacidad;                  //tamanio maximo
} Cola;

//estructura para ordenar individuos por riesgo
typedef struct IndividuoRiesgo{
    Individuo *individuo;           //puntero al individuo
    float riesgo_calculado;         //valor de riesgo calculado
} IndividuoRiesgo;

//estructura para aristas del arbol de expansion minima
typedef struct AristaMST{
    int territorio_u;               //territorio origen
    int territorio_v;               //territorio destino
    float peso;                     //peso de la arista
} AristaMST;

//estructura auxiliar para ordenamiento de individuos
typedef struct IndividuoOrden{
    Individuo *individuo;           //puntero al individuo
    float valor_orden;              //valor por el cual ordenar
} IndividuoOrden;

//enumeracion con indices de los 20 territorios
enum TerritoriosIdx {
    CHINA = 0, JAPON, KOREA, TURQUIA, ARABIA, RUSIA, GRECIA, CROACIA, HUNGRIA, POLONIA,
    FINLANDIA, SUECIA, DINAMARCA, ALEMANIA, FRANCIA, ESPANA, PORTUGAL, ITALIA, EUA, REINO_UNIDO
};

//variable global para generar ids unicos
int IDs = 0;

//declaraciones de funciones basicas
void MENU();
void CrearTerritorio(Territorio *t, int id, const char *nom, int cap);
void CrearConexiones(Mapa *grafo);
float Azar(float lim_inf, float lim_sup);
void AgregarConexion(Mapa *grafo, int t1, int t2, float peso);
void InicializarGrafo(Mapa *grafo);
void AgregarIndividuo(Territorio *territorio, Individuo *individuo);
void CrearIndividuos(Mapa *grafo, int territorio_id, sqlite3 *db);
void ConsultaSQL(int territorio_id, char *buffer);

//funciones del minheap para dijkstra y prim
MinHeap* CrearMinHeap(int capacidad);
void IntercambiarNodos(NodoHeap *a, NodoHeap *b);
void MinHeapify(MinHeap *heap, int idx);
int EstaVacio(MinHeap *heap);
NodoHeap ExtraerMin(MinHeap *heap);
void DisminuirClave(MinHeap *heap, int vertice, float nueva_prioridad);
void InsertarHeap(MinHeap *heap, int vertice, float prioridad);
void LiberarHeap(MinHeap *heap);

//funciones del sistema
void InicializarCepas(Mapa *grafo);
void GenerarRedContactos(Mapa *grafo);
void AgregarContacto(Individuo *ind1, Individuo *ind2, float prob);
int ExisteContacto(Individuo *ind, int id_otro);

//funciones de ordenamiento - o(n log n)
void OrdenarPorRiesgo(Mapa *grafo);
void OrdenarPorGrado(Mapa *grafo);
//ordena individuos por territorio de origen
void OrdenarPorTerritorio(Mapa *grafo);
//compara dos individuos por su nivel de riesgo
int CompararPorRiesgo(const void *a, const void *b);
//compara dos individuos por su grado de contactos
int CompararPorGrado(const void *a, const void *b);
//compara dos individuos por su territorio
int CompararPorTerritorio(const void *a, const void *b);
//muestra estadisticas generales del sistema
void AnalisisDatos(Mapa *grafo);

//funciones bfs para deteccion de brotes - o(v+e)
Cola* CrearCola(int capacidad);
int ColaVacia(Cola *cola);
void Encolar(Cola *cola, int item);
int Desencolar(Cola *cola);
void LiberarCola(Cola *cola);
void IniciarBrote(Mapa *grafo, int territorio_inicial, int cepa_id, int num_infectados);
void DetectarBrotes(Mapa *grafo);
void BFS_Brote(Mapa *grafo, int territorio_origen, int *visitados, int *cluster, int *tam_cluster);
void MostrarEstadisticasBrotes(Mapa *grafo);

//funciones de propagacion temporal
//simula la propagacion del virus durante varios dias
void SimularPropagacion(Mapa *grafo, int num_dias);
//avanza un dia en la simulacion procesando contagios y recuperaciones
void AvanzarUnDia(Mapa *grafo, int dia_actual);
//propaga el contagio de infectados a sus contactos sanos
void Propagar_Contagio(Mapa *grafo);
//actualiza el estado de cada individuo segun el tiempo transcurrido
void ActualizarEstados(Mapa *grafo, int dia_actual);
//muestra las estadisticas de un dia especifico
void MostrarEstadoDia(Mapa *grafo, int dia);
void GenerarReportePropagacion(Mapa *grafo, int dia_actual);
//cuenta el numero de infectados activos en todo el sistema
int ContarInfectadosActivos(Mapa *grafo);

//funciones greedy para minimizar riesgo - o(n log n)
void CalcularRiesgoIndividuos(Mapa *grafo, IndividuoRiesgo **lista_riesgo, int *total);
//selecciona individuos para vacunar usando algoritmo greedy
void MinimizarRiesgoGreedy(Mapa *grafo, int num_vacunas);
//calcula el riesgo total sumando riesgos de todos los individuos
float CalcularRiesgoTotal(Mapa *grafo);
//marca un individuo como vacunado reduciendo su riesgo a cero
void VacunarIndividuo(Individuo *ind);

//funciones prim para arbol de expansion minima - o((n+m) log n)
//encuentra el arbol de expansion minima usando prim con heap
void AlgoritmoPrim(Mapa *grafo, int territorio_inicio);
//muestra el arbol de expansion minima resultante
void MostrarMST(Mapa *grafo, int *padre, float *clave);
float CalcularPesoTotalMST(Mapa *grafo, int *padre);

//funciones hash table para busqueda en o(1)
//crea una nueva tabla hash vacia para individuos
HashTable* CrearHashTable();
int FuncionHash(int id);
//inserta un individuo en la tabla hash
void InsertarHash(HashTable *tabla, Individuo *individuo);
//busca un individuo por id en la tabla hash - o(1)
Individuo* BuscarHash(HashTable *tabla, int id);
void EliminarHash(HashTable *tabla, int id);
void MostrarEstadisticasHash(HashTable *tabla);
//inicializa la tabla hash con todos los individuos del sistema
void InicializarHashTable(Mapa *grafo);
void PruebaRendimientoHash(Mapa *grafo);
void LiberarHashTable(HashTable *tabla);

//funciones hash table para cepas - o(1)
//crea una tabla hash vacia para cepas
HashTableCepas* CrearHashTableCepas();
//calcula el indice hash para un id de cepa
int FuncionHashCepa(int id);
//inserta una cepa en la tabla hash de cepas
void InsertarHashCepa(HashTableCepas *tabla, Cepa *cepa);
//busca una cepa por id en la tabla hash - o(1)
Cepa* BuscarHashCepa(HashTableCepas *tabla, int id);
//inicializa la tabla hash con todas las cepas
void InicializarHashTableCepas(Mapa *grafo);

//funciones para semillas iniciales de contagio
//configura las 10 semillas iniciales de contagio precargadas
void InicializarSemillas(Mapa *grafo);
//aplica las semillas infectando a los individuos correspondientes
void AplicarSemillas(Mapa *grafo);
//muestra la tabla de las 10 semillas iniciales
void MostrarSemillas(Mapa *grafo);

//funciones de ordenamiento adicionales
//ordena individuos por tiempo de infeccion ascendente
void OrdenarPorTiempoInfeccion(Mapa *grafo);
//ordena individuos por nombre alfabeticamente
void OrdenarPorNombre(Mapa *grafo);
//compara dos individuos por tiempo de infeccion
int CompararPorTiempoInfeccion(const void *a, const void *b);
//compara dos individuos por su nombre alfabeticamente
int CompararPorNombre(const void *a, const void *b);

//funciones dijkstra para rutas criticas - o((n+m) log n)
//encuentra las rutas mas cortas desde un territorio usando dijkstra
void AlgoritmoDijkstra(Mapa *grafo, int territorio_origen);
void MostrarRutasCriticas(Mapa *grafo, int origen, float *distancia, int *padre);
//reconstruye y muestra la ruta desde origen hasta destino
void ReconstruirRuta(Mapa *grafo, int origen, int destino, int *padre);

//funciones trie para clustering de cepas - o(n*l)
//crea un nuevo nodo vacio para el trie
NodoTrie* CrearNodoTrie();
//crea un nuevo trie vacio
Trie* CrearTrie();
//inserta una palabra en el trie asociada a un id de cepa
void InsertarEnTrie(Trie *trie, const char *palabra, int cepa_id);
//busca una palabra en el trie y retorna el id de cepa si existe
int BuscarEnTrie(Trie *trie, const char *palabra);
//busca todas las cepas que empiezan con un prefijo dado
void BuscarPorPrefijo(Trie *trie, const char *prefijo, Mapa *grafo);
void BuscarPorPrefijoRecursivo(NodoTrie *nodo, char *prefijo_actual, int nivel, Mapa *grafo, int *resultados, int *num_resultados);
//inicializa el trie con todas las cepas del sistema
void InicializarTrie(Mapa *grafo);
void MostrarCepasPorCluster(Mapa *grafo);
void LiberarTrie(NodoTrie *nodo);
int CharAIndice(char c);

//funciones del menu
//muestra el menu de opciones de ordenamiento
void MenuOrdenamiento(Mapa *grafo);

//quicksort - o(n log n) promedio
//ordena el arreglo usando quicksort recursivo
void QuickSort(IndividuoOrden *arr, int inicio, int fin);
//particiona el arreglo alrededor de un pivote
int Particionar(IndividuoOrden *arr, int inicio, int fin);
//ordena individuos por riesgo descendente usando quicksort
void OrdenarPorRiesgoQuick(Mapa *grafo);

//quicksort para IndividuoRiesgo (usado en greedy)
void QuickSortRiesgo(IndividuoRiesgo *arr, int inicio, int fin);
int ParticionarRiesgo(IndividuoRiesgo *arr, int inicio, int fin);

//quicksort por grado (descendente)
void QuickSortGrado(IndividuoOrden *arr, int inicio, int fin);
int ParticionarGrado(IndividuoOrden *arr, int inicio, int fin);

//quicksort por territorio (ascendente)
void QuickSortTerritorio(IndividuoOrden *arr, int inicio, int fin);
int ParticionarTerritorio(IndividuoOrden *arr, int inicio, int fin);

//quicksort por tiempo de infeccion (ascendente, no infectados al final)
void QuickSortTiempoInfeccion(IndividuoOrden *arr, int inicio, int fin);
int ParticionarTiempoInfeccion(IndividuoOrden *arr, int inicio, int fin);

//quicksort por nombre (alfabetico)
void QuickSortNombre(IndividuoOrden *arr, int inicio, int fin);
int ParticionarNombre(IndividuoOrden *arr, int inicio, int fin);

//mergesort - o(n log n) garantizado
//ordena el arreglo usando mergesort recursivo
void MergeSort(IndividuoOrden *arr, int inicio, int fin);
//combina dos subarreglos ordenados en uno solo
void Merge(IndividuoOrden *arr, int inicio, int medio, int fin);
//ordena individuos por riesgo descendente usando mergesort
void OrdenarPorRiesgoMerge(Mapa *grafo);

//heapsort - o(n log n) in-place
//ordena el arreglo usando heapsort
void HeapSort(IndividuoOrden *arr, int n);
//mantiene la propiedad heap en el subarbol con raiz en i
void HeapifyOrdenamiento(IndividuoOrden *arr, int n, int i);
//ordena individuos por riesgo descendente usando heapsort
void OrdenarPorRiesgoHeap(Mapa *grafo);

//countingsort - o(n+k) lineal
//ordena por grado usando counting sort - o(n+k)
void CountingSort(Mapa *grafo);


//funcion principal del programa
int main(int argc, char const *argv[]){
    srand(42);
    
    printf("=== Iniciando BioSim ===\n");
    
    //abrir conexion a la base de datos de nombres
    sqlite3 *db;
    int rc = sqlite3_open("Nombres.db", &db);
    if(rc != SQLITE_OK){
        printf("Error al abrir base de datos\n");
        return 1;
    }
    
    Mapa mundo;
    printf("Inicializando grafo de territorios...\n");
    InicializarGrafo(&mundo);
    CrearConexiones(&mundo);
    
    printf("Creando individuos en territorios...\n");
    for(int i = 0; i < NUM_TERRITORIOS; i++){
        CrearIndividuos(&mundo, i, db);
    }
    
    printf("Inicializando %d cepas virales...\n", NUM_CEPAS);
    InicializarCepas(&mundo);
    
    printf("Generando red de contactos...\n");
    GenerarRedContactos(&mundo);
    
    printf("Inicializando 10 semillas de contagio...\n");
    InicializarSemillas(&mundo);
    AplicarSemillas(&mundo);
    
    printf("\n=== Sistema inicializado correctamente ===\n");
    printf("Territorios: %d\n", mundo.num_territorios);
    printf("Conexiones: %d\n", mundo.num_conexiones);
    
    int total_individuos = 0;
    for(int i = 0; i < NUM_TERRITORIOS; i++){
        total_individuos += mundo.territorios[i].num_individuos;
    }
    printf("Individuos totales: %d\n", total_individuos);
    printf("Cepas: %d\n", mundo.num_cepas);
    printf("Semillas iniciales: %d infectados\n", mundo.num_semillas);
    
    int opcion = -1;
    
    while(opcion != 0){
        MENU();
        printf("Seleccione una opcion: ");
        scanf("%d", &opcion);
        getchar();
        
        switch(opcion){
            case 1:
                printf("\nInicializacion y analisis de datos\n");
                
                MenuOrdenamiento(&mundo);
                
                int opcion_fase1;
                scanf("%d", &opcion_fase1);
                getchar();
                
                switch(opcion_fase1){
                    case 1:
                        OrdenarPorRiesgoQuick(&mundo);
                        break;
                    case 2:
                        OrdenarPorRiesgoMerge(&mundo);
                        break;
                    case 3:
                        OrdenarPorRiesgoHeap(&mundo);
                        break;
                    case 4:
                        CountingSort(&mundo);
                        break;
                    case 5:
                        OrdenarPorTiempoInfeccion(&mundo);
                        break;
                    case 6:
                        OrdenarPorNombre(&mundo);
                        break;
                    case 7:
                        OrdenarPorTerritorio(&mundo);
                        break;
                    case 8:
                        AnalisisDatos(&mundo);
                        break;
                    case 9:
                        MostrarSemillas(&mundo);
                        break;
                    default:
                        printf("\nOpción inválida\n");
                        break;
                }
                
                printf("\nPresione Enter para continuar...");
                getchar();
                break;
                
            case 2:
                printf("\nDeteccion de brotes con BFS\n");
                printf("\n¿Desea iniciar un nuevo brote? (1=Si, 0=No): ");
                int iniciar;
                scanf("%d", &iniciar);
                getchar();
                
                if(iniciar == 1){
                    printf("\nTerritorios disponibles:\n");
                    for(int i = 0; i < NUM_TERRITORIOS; i++){
                        printf("%2d. %s\n", i, mundo.territorios[i].Nombre);
                    }
                    
                    int terr_id, cepa_id, num_inf;
                    printf("\nTerritorio inicial (0-%d): ", NUM_TERRITORIOS-1);
                    scanf("%d", &terr_id);
                    printf("Cepa (0-%d): ", NUM_CEPAS-1);
                    scanf("%d", &cepa_id);
                    printf("Numero de infectados iniciales: ");
                    scanf("%d", &num_inf);
                    getchar();
                    
                    IniciarBrote(&mundo, terr_id, cepa_id, num_inf);
                }
                
                DetectarBrotes(&mundo);
                printf("\nPresione Enter para continuar...");
                getchar();
                break;
                
            case 3:
                printf("\nPropagacion temporal\n");
                
                int infectados_check = ContarInfectadosActivos(&mundo);
                
                if(infectados_check == 0){
                    printf("\n⚠ No hay brote activo.\n");
                    printf("Inicie un brote en Fase 2 primero.\n");
                } else {
                    printf("\nInfectados activos: %d\n", infectados_check);
                    printf("Dias a simular (recomendado: 10-30): ");
                    
                    int num_dias;
                    scanf("%d", &num_dias);
                    getchar();
                    
                    if(num_dias > 0 && num_dias <= 100){
                        SimularPropagacion(&mundo, num_dias);
                    } else {
                        printf("Numero invalido (1-100).\n");
                    }
                }
                
                printf("\nPresione Enter para continuar...");
                getchar();
                break;
                
            case 4:
                printf("\nMinimizacion del riesgo total con Greedy\n");
                
                int num_vacunas;
                printf("\nNumero de vacunas/recursos disponibles: ");
                scanf("%d", &num_vacunas);
                getchar();
                
                MinimizarRiesgoGreedy(&mundo, num_vacunas);
                
                printf("\nPresione Enter para continuar...");
                getchar();
                break;
                
            case 5:
                printf("\nIdentificacion de rutas criticas con Dijkstra\n");
                
                printf("\nTerritorios disponibles:\n");
                for(int i = 0; i < NUM_TERRITORIOS; i++){
                    printf("%2d. %s\n", i, mundo.territorios[i].Nombre);
                }
                
                int territorio_dijkstra;
                printf("\nTerritorio origen (0-%d): ", NUM_TERRITORIOS-1);
                scanf("%d", &territorio_dijkstra);
                getchar();
                
                if(territorio_dijkstra >= 0 && territorio_dijkstra < NUM_TERRITORIOS){
                    AlgoritmoDijkstra(&mundo, territorio_dijkstra);
                } else {
                    printf("Territorio invalido.\n");
                }
                
                printf("\nPresione Enter para continuar...");
                getchar();
                break;
                
            case 6:
                printf("\nCalculo de rutas optimas de contencion con Prim\n");
                
                printf("\nTerritorios disponibles:\n");
                for(int i = 0; i < NUM_TERRITORIOS; i++){
                    printf("%2d. %s\n", i, mundo.territorios[i].Nombre);
                }
                
                int territorio_prim;
                printf("\nTerritorio de inicio para MST (0-%d): ", NUM_TERRITORIOS-1);
                scanf("%d", &territorio_prim);
                getchar();
                
                AlgoritmoPrim(&mundo, territorio_prim);
                
                printf("\nPresione Enter para continuar...");
                getchar();
                break;
                
            case 7:
                printf("\nClustering de cepas similares con Trie\n");
                printf("\nOpciones:\n");
                printf("1. Inicializar Trie\n");
                printf("2. Buscar cepa por nombre exacto\n");
                printf("3. Buscar cepas por prefijo\n");
                printf("4. Mostrar clustering completo\n");
                printf("Seleccione: ");
                
                int opcion_trie;
                scanf("%d", &opcion_trie);
                getchar();
                
                switch(opcion_trie){
                    case 1:
                        InicializarTrie(&mundo);
                        break;
                        
                    case 2:
                        if(mundo.trie_cepas == NULL){
                            printf("\nPrimero debe inicializar el Trie (opcion 1)\n");
                        } else {
                            char nombre_buscar[30];
                            printf("\nNombre de la cepa a buscar: ");
                            fgets(nombre_buscar, 30, stdin);
                            nombre_buscar[strcspn(nombre_buscar, "\n")] = 0;
                            
                            int cepa_id = BuscarEnTrie(mundo.trie_cepas, nombre_buscar);
                            
                            if(cepa_id != -1){
                                printf("\n✓ CEPA ENCONTRADA:\n");
                                Cepa *cepa = &mundo.cepas[cepa_id];
                                printf("  ID: %d\n", cepa->ID);
                                printf("  Nombre: %s\n", cepa->Nombre);
                                printf("  Tasa contagio: %.3f\n", cepa->Tasa_contagio);
                                printf("  Tasa mortalidad: %.3f\n", cepa->Tasa_mortalidad);
                                printf("  Tiempo incubación: %d días\n", cepa->Tiempo_incubacion);
                                printf("  Tiempo recuperación: %d días\n", cepa->Tiempo_recuperacion);
                            } else {
                                printf("\n✗ Cepa no encontrada\n");
                            }
                        }
                        break;
                        
                    case 3:
                        if(mundo.trie_cepas == NULL){
                            printf("\nPrimero debe inicializar el Trie (opcion 1)\n");
                        } else {
                            char prefijo[30];
                            printf("\nPrefijo a buscar (ej: Alpha, Beta, Gamma): ");
                            fgets(prefijo, 30, stdin);
                            prefijo[strcspn(prefijo, "\n")] = 0;
                            
                            BuscarPorPrefijo(mundo.trie_cepas, prefijo, &mundo);
                        }
                        break;
                        
                    case 4:
                        MostrarCepasPorCluster(&mundo);
                        break;
                        
                    default:
                        printf("\nOpcion invalida\n");
                        break;
                }
                
                printf("\nPresione Enter para continuar...");
                getchar();
                break;
                
            case 8:
                printf("\nAlmacenamiento eficiente y consulta rapida O(1)\n");
                printf("\nOpciones:\n");
                printf("--- Hash Table Individuos ---\n");
                printf("1. Inicializar Hash Table Individuos\n");
                printf("2. Buscar individuo por ID - O(1)\n");
                printf("3. Mostrar estadisticas Hash Individuos\n");
                printf("4. Prueba de rendimiento\n");
                printf("--- Hash Table Cepas ---\n");
                printf("5. Inicializar Hash Table Cepas\n");
                printf("6. Buscar cepa por ID - O(1)\n");
                printf("Seleccione: ");
                
                int opcion_hash;
                scanf("%d", &opcion_hash);
                getchar();
                
                switch(opcion_hash){
                    case 1:
                        InicializarHashTable(&mundo);
                        break;
                        
                    case 2:
                        if(mundo.hash_individuos == NULL){
                            printf("\nPrimero debe inicializar la Hash Table (opcion 1)\n");
                        } else {
                            int id_buscar;
                            printf("\nID del individuo a buscar: ");
                            scanf("%d", &id_buscar);
                            getchar();
                            
                            Individuo *resultado = BuscarHash(mundo.hash_individuos, id_buscar);
                            
                            if(resultado != NULL){
                                printf("\n✓ INDIVIDUO ENCONTRADO en O(1):\n");
                                printf("  ID: %d\n", resultado->ID);
                                printf("  Nombre: %s\n", resultado->Nombre);
                                printf("  Territorio: %s\n", mundo.territorios[resultado->Territorio_ID].Nombre);
                                printf("  Riesgo inicial: %.3f\n", resultado->Riesgo_inicial);
                                printf("  Tiempo infección: %d\n", resultado->t_infeccion);
                                printf("  Estado: %s\n", 
                                    resultado->Infectado ? "INFECTADO" : 
                                    resultado->Recuperado ? "RECUPERADO" :
                                    resultado->Fallecido ? "FALLECIDO" : "SANO");
                                
                                if(resultado->Infectado && resultado->Cepa_ID >= 0){
                                    printf("  Cepa: %s\n", mundo.cepas[resultado->Cepa_ID].Nombre);
                                }
                                
                                int num_contactos = 0;
                                Contacto *c = resultado->contactos;
                                while(c != NULL){
                                    num_contactos++;
                                    c = c->sgt;
                                }
                                printf("  Numero de contactos: %d\n", num_contactos);
                            } else {
                                printf("\n✗ Individuo con ID %d no encontrado\n", id_buscar);
                            }
                        }
                        break;
                        
                    case 3:
                        if(mundo.hash_individuos == NULL){
                            printf("\nPrimero debe inicializar la Hash Table (opcion 1)\n");
                        } else {
                            MostrarEstadisticasHash(mundo.hash_individuos);
                        }
                        break;
                        
                    case 4:
                        if(mundo.hash_individuos == NULL){
                            printf("\nPrimero debe inicializar la Hash Table (opcion 1)\n");
                        } else {
                            PruebaRendimientoHash(&mundo);
                        }
                        break;
                    
                    case 5:
                        InicializarHashTableCepas(&mundo);
                        break;
                        
                    case 6:
                        if(mundo.hash_cepas == NULL){
                            printf("\nPrimero debe inicializar la Hash Table de Cepas (opcion 5)\n");
                        } else {
                            int cepa_id;
                            printf("\nID de la cepa a buscar (0-%d): ", NUM_CEPAS-1);
                            scanf("%d", &cepa_id);
                            getchar();
                            
                            Cepa *cepa_encontrada = BuscarHashCepa(mundo.hash_cepas, cepa_id);
                            
                            if(cepa_encontrada != NULL){
                                printf("\n✓ CEPA ENCONTRADA en O(1):\n");
                                printf("  ID: %d\n", cepa_encontrada->ID);
                                printf("  Nombre: %s\n", cepa_encontrada->Nombre);
                                printf("  Tasa contagio (beta): %.3f\n", cepa_encontrada->Tasa_contagio);
                                printf("  Letalidad: %.3f\n", cepa_encontrada->Tasa_mortalidad);
                                printf("  Tiempo incubación: %d días\n", cepa_encontrada->Tiempo_incubacion);
                                printf("  Tiempo recuperación (gamma): %d días\n", cepa_encontrada->Tiempo_recuperacion);
                            } else {
                                printf("\n✗ Cepa con ID %d no encontrada\n", cepa_id);
                            }
                        }
                        break;
                        
                    default:
                        printf("\nOpcion invalida\n");
                        break;
                }
                
                printf("\nPresione Enter para continuar...");
                getchar();
                break;
                
            case 0:
                printf("\nSaliendo de BioSim...\n");
                break;
                
            default:
                printf("\nOpcion invalida.\n");
                printf("Presione Enter para continuar...");
                getchar();
                break;
        }
    }
    
    sqlite3_close(db);
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
    grafo->hash_individuos = NULL;
    grafo->hash_cepas = NULL;
    grafo->trie_cepas = NULL;
    grafo->num_semillas = 0;
    
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

//muestra el menu principal con las 8 fases del proyecto
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

//agrega un individuo a un territorio si hay espacio disponible
void AgregarIndividuo(Territorio *territorio, Individuo *individuo){
    if(territorio->num_individuos < MAX_INDIVIDUOS){
        territorio->individuos[territorio->num_individuos] = individuo;
        territorio->num_individuos++;
        individuo->Territorio_ID = territorio->ID;
    }
}

//crea los individuos de un territorio leyendo nombres de la base de datos
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
    
    //mapear territorio_id a nombre de pais en la base de datos
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

//genera un numero aleatorio entre lim_inf y lim_sup
float Azar(float lim_inf, float lim_sup){ 
    float num = rand() / (float) RAND_MAX;
    num = lim_inf + num * (lim_sup - lim_inf);
    return num;
}

//=============================================================
//funciones del minheap para dijkstra y prim
//=============================================================

//crea un heap minimo con la capacidad especificada
MinHeap* CrearMinHeap(int capacidad){
    MinHeap *heap = (MinHeap*)malloc(sizeof(MinHeap));
    heap->capacidad = capacidad;
    heap->size = 0;
    heap->nodos = (NodoHeap*)malloc(capacidad * sizeof(NodoHeap));
    heap->posiciones = (int*)malloc(capacidad * sizeof(int));
    
    for(int i = 0; i < capacidad; i++){
        heap->posiciones[i] = -1;
    }
    
    return heap;
}

//intercambia dos nodos del heap
void IntercambiarNodos(NodoHeap *a, NodoHeap *b){
    NodoHeap temp = *a;
    *a = *b;
    *b = temp;
}

//mantiene la propiedad de heap minimo bajando el nodo en idx
void MinHeapify(MinHeap *heap, int idx){
    int smallest = idx;
    int left = 2 * idx + 1;
    int right = 2 * idx + 2;
    
    if(left < heap->size && heap->nodos[left].prioridad < heap->nodos[smallest].prioridad)
        smallest = left;
    
    if(right < heap->size && heap->nodos[right].prioridad < heap->nodos[smallest].prioridad)
        smallest = right;
    
    if(smallest != idx){
        heap->posiciones[heap->nodos[smallest].vertice] = idx;
        heap->posiciones[heap->nodos[idx].vertice] = smallest;
        
        IntercambiarNodos(&heap->nodos[smallest], &heap->nodos[idx]);
        MinHeapify(heap, smallest);
    }
}

//verifica si el heap esta vacio
int EstaVacio(MinHeap *heap){
    return heap->size == 0;
}

NodoHeap ExtraerMin(MinHeap *heap){
    if(EstaVacio(heap)){
        NodoHeap nodo_vacio = {-1, -1.0};
        return nodo_vacio;
    }
    
    NodoHeap raiz = heap->nodos[0];
    NodoHeap ultimo = heap->nodos[heap->size - 1];
    heap->nodos[0] = ultimo;
    
    heap->posiciones[raiz.vertice] = -1;
    heap->posiciones[ultimo.vertice] = 0;
    
    heap->size--;
    MinHeapify(heap, 0);
    
    return raiz;
}

//disminuye la prioridad de un vertice y lo sube en el heap
void DisminuirClave(MinHeap *heap, int vertice, float nueva_prioridad){
    int i = heap->posiciones[vertice];
    
    if(i == -1) return;
    
    heap->nodos[i].prioridad = nueva_prioridad;
    
    while(i > 0 && heap->nodos[i].prioridad < heap->nodos[(i-1)/2].prioridad){
        heap->posiciones[heap->nodos[i].vertice] = (i-1)/2;
        heap->posiciones[heap->nodos[(i-1)/2].vertice] = i;
        
        IntercambiarNodos(&heap->nodos[i], &heap->nodos[(i-1)/2]);
        i = (i-1)/2;
    }
}

//inserta un nuevo vertice con su prioridad en el heap
void InsertarHeap(MinHeap *heap, int vertice, float prioridad){
    if(heap->size == heap->capacidad) return;
    
    heap->size++;
    int i = heap->size - 1;
    
    heap->nodos[i].vertice = vertice;
    heap->nodos[i].prioridad = prioridad;
    heap->posiciones[vertice] = i;
    
    while(i > 0 && heap->nodos[i].prioridad < heap->nodos[(i-1)/2].prioridad){
        heap->posiciones[heap->nodos[i].vertice] = (i-1)/2;
        heap->posiciones[heap->nodos[(i-1)/2].vertice] = i;
        
        IntercambiarNodos(&heap->nodos[i], &heap->nodos[(i-1)/2]);
        i = (i-1)/2;
    }
}

//libera la memoria del heap
void LiberarHeap(MinHeap *heap){
    free(heap->nodos);
    free(heap->posiciones);
    free(heap);
}

//=============================================================
//inicializacion de las 50 cepas virales
//=============================================================

//crea las 50 cepas con nombres y caracteristicas aleatorias
void InicializarCepas(Mapa *grafo){
    const char *prefijos[] = {"Alpha", "Beta", "Gamma", "Delta", "Epsilon", 
                              "Zeta", "Eta", "Theta", "Iota", "Kappa"};
    const char *sufijos[] = {"Flu", "Pox", "Fever", "Virus", "Strain"};
    
    grafo->num_cepas = NUM_CEPAS;
    
    for(int i = 0; i < NUM_CEPAS; i++){
        grafo->cepas[i].ID = i;
        
        sprintf(grafo->cepas[i].Nombre, "%s-%s-%d", 
                prefijos[i % 10], 
                sufijos[i % 5], 
                i + 1);
        
        grafo->cepas[i].Tasa_contagio = Azar(0.1, 0.9);
        grafo->cepas[i].Tasa_mortalidad = Azar(0.01, 0.15);
        grafo->cepas[i].Tiempo_incubacion = (int)Azar(1, 14);
        grafo->cepas[i].Tiempo_recuperacion = (int)Azar(7, 21);
        
        strcpy(grafo->cepas[i].Sintomas, "Fiebre, tos, fatiga");
    }
}

//=============================================================
//generacion de la red de contactos entre individuos
//=============================================================

//agrega un contacto bidireccional entre dos individuos
void AgregarContacto(Individuo *ind1, Individuo *ind2, float prob){
    if(ind1 == NULL || ind2 == NULL) return;
    
    Contacto *nuevo = (Contacto*)malloc(sizeof(Contacto));
    nuevo->u_individuo = ind1->ID;
    nuevo->v_individuo = ind2->ID;
    nuevo->prob_contagio = prob;
    nuevo->sgt = ind1->contactos;
    ind1->contactos = nuevo;
    
    Contacto *nuevo2 = (Contacto*)malloc(sizeof(Contacto));
    nuevo2->u_individuo = ind2->ID;
    nuevo2->v_individuo = ind1->ID;
    nuevo2->prob_contagio = prob;
    nuevo2->sgt = ind2->contactos;
    ind2->contactos = nuevo2;
}

//verifica si ya existe un contacto entre dos individuos
int ExisteContacto(Individuo *ind, int id_otro){
    if(ind == NULL) return 0;
    
    Contacto *actual = ind->contactos;
    while(actual != NULL){
        if(actual->v_individuo == id_otro) return 1;
        actual = actual->sgt;
    }
    return 0;
}

void GenerarRedContactos(Mapa *grafo){
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        
        if(territorio->num_individuos == 0) continue;
        
        for(int i = 0; i < territorio->num_individuos; i++){
            Individuo *ind1 = territorio->individuos[i];
            if(ind1 == NULL) continue;
            
            int num_contactos = (int)Azar(2, ind1->Grado_inicial);
            if(num_contactos > territorio->num_individuos - 1){
                num_contactos = territorio->num_individuos - 1;
            }
            
            int contactos_creados = 0;
            
            for(int intentos = 0; intentos < num_contactos * 3 && contactos_creados < num_contactos; intentos++){
                int idx_aleatorio = (int)Azar(0, territorio->num_individuos);
                
                if(idx_aleatorio >= territorio->num_individuos) continue;
                if(idx_aleatorio == i) continue;
                
                Individuo *ind2 = territorio->individuos[idx_aleatorio];
                if(ind2 == NULL) continue;
                
                if(!ExisteContacto(ind1, ind2->ID)){
                    float prob = Azar(0.1, 0.8);
                    AgregarContacto(ind1, ind2, prob);
                    contactos_creados++;
                }
            }
        }
        
        for(int t2 = t + 1; t2 < NUM_TERRITORIOS; t2++){
            if(grafo->matrix[t][t2] > 0.0){
                Territorio *territorio2 = &grafo->territorios[t2];
                
                if(territorio2->num_individuos == 0) continue;
                
                int contactos_inter = (int)Azar(1, 3);
                
                for(int c = 0; c < contactos_inter; c++){
                    int idx1 = (int)Azar(0, territorio->num_individuos);
                    int idx2 = (int)Azar(0, territorio2->num_individuos);
                    
                    if(idx1 >= territorio->num_individuos || idx2 >= territorio2->num_individuos) continue;
                    
                    Individuo *ind1 = territorio->individuos[idx1];
                    Individuo *ind2 = territorio2->individuos[idx2];
                    
                    if(ind1 == NULL || ind2 == NULL) continue;
                    
                    if(!ExisteContacto(ind1, ind2->ID)){
                        float prob = Azar(0.05, 0.3) * grafo->matrix[t][t2];
                        AgregarContacto(ind1, ind2, prob);
                    }
                }
            }
        }
    }
}

//=============================================================
//deteccion de brotes usando bfs - o(v+e)
//=============================================================

//crea una cola con la capacidad especificada
Cola* CrearCola(int capacidad){
    Cola *cola = (Cola*)malloc(sizeof(Cola));
    cola->capacidad = capacidad;
    cola->frente = 0;
    cola->final = -1;
    cola->items = (int*)malloc(capacidad * sizeof(int));
    return cola;
}

//verifica si la cola esta vacia
int ColaVacia(Cola *cola){
    return cola->final < cola->frente;
}

//agrega un elemento al final de la cola
void Encolar(Cola *cola, int item){
    if(cola->final < cola->capacidad - 1){
        cola->final++;
        cola->items[cola->final] = item;
    }
}

//remueve y retorna el elemento del frente de la cola
int Desencolar(Cola *cola){
    if(ColaVacia(cola)) return -1;
    int item = cola->items[cola->frente];
    cola->frente++;
    return item;
}

//libera la memoria de la cola
void LiberarCola(Cola *cola){
    free(cola->items);
    free(cola);
}

//inicia un brote infectando individuos en un territorio
void IniciarBrote(Mapa *grafo, int territorio_inicial, int cepa_id, int num_infectados){
    if(territorio_inicial < 0 || territorio_inicial >= NUM_TERRITORIOS) return;
    if(cepa_id < 0 || cepa_id >= NUM_CEPAS) return;
    
    Territorio *territorio = &grafo->territorios[territorio_inicial];
    
    if(num_infectados > territorio->num_individuos){
        num_infectados = territorio->num_individuos;
    }
    
    printf("\n>>> Iniciando brote en %s con cepa %s <<<\n", 
           territorio->Nombre, grafo->cepas[cepa_id].Nombre);
    
    int infectados = 0;
    for(int i = 0; i < territorio->num_individuos && infectados < num_infectados; i++){
        Individuo *ind = territorio->individuos[i];
        if(ind != NULL && !ind->Infectado){
            ind->Infectado = 1;
            ind->t_infeccion = 0;
            ind->Cepa_ID = cepa_id;
            infectados++;
        }
    }
    
    printf(">>> %d infectados <<<\n", infectados);
}

//realiza bfs desde un territorio para encontrar territorios conectados con infectados
void BFS_Brote(Mapa *grafo, int territorio_origen, int *visitados, int *cluster, int *tam_cluster){
    Cola *cola = CrearCola(NUM_TERRITORIOS);
    
    visitados[territorio_origen] = 1;
    cluster[*tam_cluster] = territorio_origen;
    (*tam_cluster)++;
    
    Encolar(cola, territorio_origen);
    
    while(!ColaVacia(cola)){
        int t_actual = Desencolar(cola);
        
        for(int t = 0; t < NUM_TERRITORIOS; t++){
            if(!visitados[t] && grafo->matrix[t_actual][t] > 0.0){
                Territorio *terr = &grafo->territorios[t];
                int tiene_infectados = 0;
                
                for(int i = 0; i < terr->num_individuos; i++){
                    if(terr->individuos[i] != NULL && terr->individuos[i]->Infectado){
                        tiene_infectados = 1;
                        break;
                    }
                }
                
                if(tiene_infectados){
                    visitados[t] = 1;
                    cluster[*tam_cluster] = t;
                    (*tam_cluster)++;
                    Encolar(cola, t);
                }
            }
        }
    }
    
    LiberarCola(cola);
}

//detecta y muestra todos los brotes activos usando bfs
void DetectarBrotes(Mapa *grafo){
    printf("\n========== DETECCIÓN BFS ==========\n");
    
    int visitados[NUM_TERRITORIOS] = {0};
    int num_clusters = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        if(!visitados[t]){
            Territorio *territorio = &grafo->territorios[t];
            int tiene_infectados = 0;
            
            for(int i = 0; i < territorio->num_individuos; i++){
                if(territorio->individuos[i] != NULL && territorio->individuos[i]->Infectado){
                    tiene_infectados = 1;
                    break;
                }
            }
            
            if(tiene_infectados){
                int cluster[NUM_TERRITORIOS];
                int tam_cluster = 0;
                
                BFS_Brote(grafo, t, visitados, cluster, &tam_cluster);
                
                num_clusters++;
                printf("\n--- CLUSTER #%d ---\n", num_clusters);
                printf("Territorios (%d):\n", tam_cluster);
                
                for(int i = 0; i < tam_cluster; i++){
                    Territorio *terr = &grafo->territorios[cluster[i]];
                    
                    int infectados = 0;
                    for(int j = 0; j < terr->num_individuos; j++){
                        if(terr->individuos[j] != NULL && terr->individuos[j]->Infectado){
                            infectados++;
                        }
                    }
                    
                    printf("  - %s: %d infectados\n", terr->Nombre, infectados);
                }
            }
        }
    }
    
    if(num_clusters == 0){
        printf("\nNo hay brotes.\n");
    } else {
        printf("\n===================================\n");
        printf("TOTAL CLUSTERS: %d\n", num_clusters);
        printf("===================================\n");
    }
}

void MostrarEstadisticasBrotes(Mapa *grafo){
    int total_infectados = 0;
    int total_recuperados = 0;
    int total_fallecidos = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            Individuo *ind = territorio->individuos[i];
            if(ind != NULL){
                if(ind->Infectado) total_infectados++;
                if(ind->Recuperado) total_recuperados++;
                if(ind->Fallecido) total_fallecidos++;
            }
        }
    }
    
    printf("\n========== ESTADÍSTICAS ==========\n");
    printf("Infectados: %d\n", total_infectados);
    printf("Recuperados: %d\n", total_recuperados);
    printf("Fallecidos: %d\n", total_fallecidos);
    printf("==================================\n");
}

//=============================================================
//propagacion temporal del contagio
//=============================================================

//cuenta el numero de infectados activos en todo el sistema
int ContarInfectadosActivos(Mapa *grafo){
    int infectados = 0;
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL && territorio->individuos[i]->Infectado){
                infectados++;
            }
        }
    }
    return infectados;
}

//propaga el contagio de infectados a sus contactos sanos
void Propagar_Contagio(Mapa *grafo){
    #define MAX_CONTAGIOS_DIA 100
    
    int num_nuevos = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS && num_nuevos < MAX_CONTAGIOS_DIA; t++){
        Territorio *territorio = &grafo->territorios[t];
        
        for(int i = 0; i < territorio->num_individuos && num_nuevos < MAX_CONTAGIOS_DIA; i++){
            Individuo *infectado = territorio->individuos[i];
            
            if(infectado == NULL || !infectado->Infectado || infectado->Fallecido) continue;
            
            Contacto *c = infectado->contactos;
            int contactos_revisados = 0;
            
            while(c != NULL && contactos_revisados < 5 && num_nuevos < MAX_CONTAGIOS_DIA){
                contactos_revisados++;
                
                Individuo *contacto = NULL;
                
                for(int j = 0; j < territorio->num_individuos; j++){
                    if(territorio->individuos[j] != NULL && 
                       territorio->individuos[j]->ID == c->v_individuo){
                        contacto = territorio->individuos[j];
                        break;
                    }
                }
                
                if(contacto == NULL){
                    for(int t2 = 0; t2 < 3 && contacto == NULL; t2++){
                        if(t2 == t) continue;
                        Territorio *terr2 = &grafo->territorios[t2];
                        for(int j = 0; j < terr2->num_individuos; j++){
                            if(terr2->individuos[j] != NULL && 
                               terr2->individuos[j]->ID == c->v_individuo){
                                contacto = terr2->individuos[j];
                                break;
                            }
                        }
                    }
                }
                
                if(contacto != NULL && !contacto->Infectado && 
                   !contacto->Recuperado && !contacto->Fallecido){
                    
                    float prob_contagio = c->prob_contagio * 0.15;
                    
                    if(infectado->Cepa_ID >= 0 && infectado->Cepa_ID < NUM_CEPAS){
                        prob_contagio *= grafo->cepas[infectado->Cepa_ID].Tasa_contagio;
                    }
                    
                    if(Azar(0.0, 1.0) < prob_contagio){
                        contacto->Infectado = 1;
                        contacto->t_infeccion = 0;
                        contacto->Cepa_ID = infectado->Cepa_ID;
                        num_nuevos++;
                    }
                }
                
                c = c->sgt;
            }
        }
    }
}

//actualiza el estado de cada individuo segun el tiempo transcurrido
void ActualizarEstados(Mapa *grafo, int dia_actual){
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            Individuo *ind = territorio->individuos[i];
            
            if(ind != NULL && ind->Infectado && !ind->Fallecido){
                ind->t_infeccion++;
                
                if(ind->Cepa_ID >= 0 && ind->Cepa_ID < NUM_CEPAS){
                    Cepa *cepa = &grafo->cepas[ind->Cepa_ID];
                    
                    if(ind->t_infeccion >= cepa->Tiempo_recuperacion){
                        ind->Infectado = 0;
                        ind->Recuperado = 1;
                        ind->Riesgo_inicial = 0.0;
                    }
                    else if(ind->t_infeccion > cepa->Tiempo_incubacion){
                        float prob_muerte = cepa->Tasa_mortalidad * 0.02;
                        if(Azar(0.0, 1.0) < prob_muerte){
                            ind->Infectado = 0;
                            ind->Fallecido = 1;
                        }
                    }
                }
            }
        }
    }
}

//muestra las estadisticas de un dia especifico
void MostrarEstadoDia(Mapa *grafo, int dia){
    int total_infectados = 0;
    int total_recuperados = 0;
    int total_fallecidos = 0;
    int total_sanos = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            Individuo *ind = territorio->individuos[i];
            if(ind != NULL){
                if(ind->Infectado) total_infectados++;
                else if(ind->Recuperado) total_recuperados++;
                else if(ind->Fallecido) total_fallecidos++;
                else total_sanos++;
            }
        }
    }
    
    printf("\n[Día %3d] Sanos: %4d | Infectados: %4d | Recuperados: %4d | Fallecidos: %4d\n", 
           dia, total_sanos, total_infectados, total_recuperados, total_fallecidos);
}

void GenerarReportePropagacion(Mapa *grafo, int dia_actual){
    printf("\n========== REPORTE FINAL ==========\n");
    
    printf("Territorios afectados:\n");
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        int infectados = 0;
        int recuperados = 0;
        
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            Individuo *ind = territorio->individuos[i];
            if(ind != NULL){
                if(ind->Infectado) infectados++;
                if(ind->Recuperado) recuperados++;
            }
        }
        
        if(infectados > 0 || recuperados > 0){
            printf("%-20s: I=%d R=%d\n", territorio->Nombre, infectados, recuperados);
        }
    }
    
    printf("===================================\n");
}

//avanza un dia en la simulacion procesando contagios y recuperaciones
void AvanzarUnDia(Mapa *grafo, int dia_actual){
    Propagar_Contagio(grafo);
    ActualizarEstados(grafo, dia_actual);
}

//simula la propagacion del virus durante varios dias
void SimularPropagacion(Mapa *grafo, int num_dias){
    printf("\n========== SIMULACIÓN ==========\n");
    
    int infectados_iniciales = ContarInfectadosActivos(grafo);
    
    if(infectados_iniciales == 0){
        printf("\nNo hay brote activo.\n");
        return;
    }
    
    printf("Infectados iniciales: %d\n", infectados_iniciales);
    printf("Días: %d\n\n", num_dias);
    
    MostrarEstadoDia(grafo, 0);
    
    for(int dia = 1; dia <= num_dias; dia++){
        AvanzarUnDia(grafo, dia);
        
        if(dia % 2 == 0 || dia == num_dias){
            MostrarEstadoDia(grafo, dia);
        }
        
        int infectados = ContarInfectadosActivos(grafo);
        if(infectados == 0){
            printf("\n✓ Epidemia extinguida (día %d)\n", dia);
            MostrarEstadoDia(grafo, dia);
            break;
        }
    }
    
    GenerarReportePropagacion(grafo, num_dias);
}

//=============================================================
//minimizacion de riesgo con algoritmo greedy - o(n log n)
//=============================================================

//calcula el riesgo total sumando riesgos de todos los individuos
float CalcularRiesgoTotal(Mapa *grafo){
    float riesgo_total = 0.0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            Individuo *ind = territorio->individuos[i];
            if(ind != NULL && !ind->Infectado && !ind->Recuperado){
                int num_contactos = 0;
                Contacto *c = ind->contactos;
                while(c != NULL){
                    num_contactos++;
                    c = c->sgt;
                }
                
                float riesgo = ind->Riesgo_inicial * (1.0 + num_contactos * 0.1);
                riesgo_total += riesgo;
            }
        }
    }
    
    return riesgo_total;
}

void CalcularRiesgoIndividuos(Mapa *grafo, IndividuoRiesgo **lista_riesgo, int *total){
    *total = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            Individuo *ind = territorio->individuos[i];
            if(ind != NULL && !ind->Infectado && !ind->Recuperado && !ind->Fallecido){
                (*total)++;
            }
        }
    }
    
    *lista_riesgo = (IndividuoRiesgo*)malloc((*total) * sizeof(IndividuoRiesgo));
    int idx = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            Individuo *ind = territorio->individuos[i];
            if(ind != NULL && !ind->Infectado && !ind->Recuperado && !ind->Fallecido){
                int num_contactos = 0;
                Contacto *c = ind->contactos;
                while(c != NULL){
                    num_contactos++;
                    c = c->sgt;
                }
                
                float riesgo = ind->Riesgo_inicial * 0.3;
                riesgo += (num_contactos * 0.05);
                
                (*lista_riesgo)[idx].individuo = ind;
                (*lista_riesgo)[idx].riesgo_calculado = riesgo;
                idx++;
            }
        }
    }
}

int CompararRiesgo(const void *a, const void *b){
    IndividuoRiesgo *ir_a = (IndividuoRiesgo*)a;
    IndividuoRiesgo *ir_b = (IndividuoRiesgo*)b;
    
    if(ir_b->riesgo_calculado > ir_a->riesgo_calculado) return 1;
    if(ir_b->riesgo_calculado < ir_a->riesgo_calculado) return -1;
    return 0;
}

//marca un individuo como vacunado reduciendo su riesgo a cero
void VacunarIndividuo(Individuo *ind){
    ind->Recuperado = 1;
    ind->Riesgo_inicial = 0.0;
}

//selecciona individuos para vacunar usando algoritmo greedy
void MinimizarRiesgoGreedy(Mapa *grafo, int num_vacunas){
    printf("\n========== GREEDY ==========\n");
    
    float riesgo_inicial = CalcularRiesgoTotal(grafo);
    printf("Riesgo inicial: %.2f\n", riesgo_inicial);
    
    IndividuoRiesgo *lista_riesgo;
    int total_individuos;
    
    CalcularRiesgoIndividuos(grafo, &lista_riesgo, &total_individuos);
    
    if(num_vacunas > total_individuos){
        num_vacunas = total_individuos;
    }
    
    printf("Total sanos: %d\n", total_individuos);
    printf("Vacunas: %d\n", num_vacunas);
    
    //ordenar por riesgo usando quicksort propio
    QuickSortRiesgo(lista_riesgo, 0, total_individuos - 1);
    
    printf("\nTop 10 mayor riesgo:\n");
    for(int i = 0; i < 10 && i < total_individuos; i++){
        printf("%2d. %s (%s) - %.3f\n", 
               i+1,
               lista_riesgo[i].individuo->Nombre,
               grafo->territorios[lista_riesgo[i].individuo->Territorio_ID].Nombre,
               lista_riesgo[i].riesgo_calculado);
    }
    
    printf("\nVacunando %d individuos...\n", num_vacunas);
    for(int i = 0; i < num_vacunas; i++){
        VacunarIndividuo(lista_riesgo[i].individuo);
    }
    
    float riesgo_final = CalcularRiesgoTotal(grafo);
    float reduccion = ((riesgo_inicial - riesgo_final) / riesgo_inicial) * 100.0;
    
    printf("\n===============================\n");
    printf("Riesgo inicial: %.2f\n", riesgo_inicial);
    printf("Riesgo final:   %.2f\n", riesgo_final);
    printf("Reducción:      %.2f%%\n", reduccion);
    printf("===============================\n");
    
    free(lista_riesgo);
}

//=============================================================
//arbol de expansion minima con prim - o((n+m) log n)
//=============================================================

float CalcularPesoTotalMST(Mapa *grafo, int *padre){
    float peso_total = 0.0;
    
    for(int v = 1; v < NUM_TERRITORIOS; v++){
        if(padre[v] != -1){
            peso_total += grafo->matrix[padre[v]][v];
        }
    }
    
    return peso_total;
}

//muestra el arbol de expansion minima resultante
void MostrarMST(Mapa *grafo, int *padre, float *clave){
    printf("\n========== MST ==========\n");
    printf("Conexiones:\n\n");
    
    for(int v = 1; v < NUM_TERRITORIOS; v++){
        if(padre[v] != -1){
            printf("%-20s <-> %-20s (%.3f)\n",
                   grafo->territorios[padre[v]].Nombre,
                   grafo->territorios[v].Nombre,
                   grafo->matrix[padre[v]][v]);
        }
    }
    
    float peso_total = CalcularPesoTotalMST(grafo, padre);
    
    printf("\n=========================\n");
    printf("PESO TOTAL: %.3f\n", peso_total);
    printf("ARISTAS: %d\n", NUM_TERRITORIOS - 1);
    printf("=========================\n");
}

//encuentra el arbol de expansion minima usando prim con heap
void AlgoritmoPrim(Mapa *grafo, int territorio_inicio){
    printf("\n========== PRIM ==========\n");
    printf("Desde: %s\n", grafo->territorios[territorio_inicio].Nombre);
    
    if(territorio_inicio < 0 || territorio_inicio >= NUM_TERRITORIOS){
        printf("Territorio inválido\n");
        return;
    }
    
    int *padre = (int*)malloc(NUM_TERRITORIOS * sizeof(int));
    float *clave = (float*)malloc(NUM_TERRITORIOS * sizeof(float));
    int *en_mst = (int*)malloc(NUM_TERRITORIOS * sizeof(int));
    
    for(int i = 0; i < NUM_TERRITORIOS; i++){
        padre[i] = -1;
        clave[i] = 999999.0;
        en_mst[i] = 0;
    }
    
    clave[territorio_inicio] = 0.0;
    
    MinHeap *heap = CrearMinHeap(NUM_TERRITORIOS);
    
    for(int v = 0; v < NUM_TERRITORIOS; v++){
        InsertarHeap(heap, v, clave[v]);
    }
    
    int vertices_procesados = 0;
    
    while(!EstaVacio(heap) && vertices_procesados < NUM_TERRITORIOS){
        NodoHeap min_nodo = ExtraerMin(heap);
        int u = min_nodo.vertice;
        
        if(u == -1) break;
        
        en_mst[u] = 1;
        vertices_procesados++;
        
        for(int v = 0; v < NUM_TERRITORIOS; v++){
            if(grafo->matrix[u][v] > 0.0 && !en_mst[v]){
                float prioridad = 1.0 - grafo->matrix[u][v];
                
                if(prioridad < clave[v]){
                    padre[v] = u;
                    clave[v] = prioridad;
                    DisminuirClave(heap, v, prioridad);
                }
            }
        }
    }
    
    MostrarMST(grafo, padre, clave);
    
    LiberarHeap(heap);
    free(padre);
    free(clave);
    free(en_mst);
}

//=============================================================
//tabla hash para busqueda de individuos en o(1)
//=============================================================

//crea una nueva tabla hash vacia para individuos
HashTable* CrearHashTable(){
    HashTable *tabla = (HashTable*)malloc(sizeof(HashTable));
    tabla->num_elementos = 0;
    
    for(int i = 0; i < HASH_SIZE; i++){
        tabla->tabla[i] = NULL;
    }
    
    return tabla;
}

//calcula el indice hash para un id dado
int FuncionHash(int id){
    return id % HASH_SIZE;
}

//inserta un individuo en la tabla hash
void InsertarHash(HashTable *tabla, Individuo *individuo){
    if(individuo == NULL) return;
    
    int indice = FuncionHash(individuo->ID);
    
    NodoHash *nuevo = (NodoHash*)malloc(sizeof(NodoHash));
    nuevo->ID = individuo->ID;
    nuevo->individuo = individuo;
    nuevo->siguiente = tabla->tabla[indice];
    
    tabla->tabla[indice] = nuevo;
    tabla->num_elementos++;
}

//busca un individuo por id en la tabla hash - o(1)
Individuo* BuscarHash(HashTable *tabla, int id){
    int indice = FuncionHash(id);
    
    NodoHash *actual = tabla->tabla[indice];
    
    while(actual != NULL){
        if(actual->ID == id){
            return actual->individuo;
        }
        actual = actual->siguiente;
    }
    
    return NULL;
}

void EliminarHash(HashTable *tabla, int id){
    int indice = FuncionHash(id);
    
    NodoHash *actual = tabla->tabla[indice];
    NodoHash *anterior = NULL;
    
    while(actual != NULL){
        if(actual->ID == id){
            if(anterior == NULL){
                tabla->tabla[indice] = actual->siguiente;
            } else {
                anterior->siguiente = actual->siguiente;
            }
            
            free(actual);
            tabla->num_elementos--;
            return;
        }
        
        anterior = actual;
        actual = actual->siguiente;
    }
}

void MostrarEstadisticasHash(HashTable *tabla){
    printf("\n========== HASH TABLE ==========\n");
    printf("Tamaño: %d\n", HASH_SIZE);
    printf("Elementos: %d\n", tabla->num_elementos);
    
    if(tabla->num_elementos == 0){
        printf("\n⚠ Tabla vacía.\n");
        printf("================================\n");
        return;
    }
    
    printf("Factor de carga: %.4f\n", (float)tabla->num_elementos / HASH_SIZE);
    
    int buckets_usados = 0;
    int max_colisiones = 0;
    int total_colisiones = 0;
    
    for(int i = 0; i < HASH_SIZE; i++){
        if(tabla->tabla[i] != NULL){
            buckets_usados++;
            
            int longitud = 0;
            NodoHash *actual = tabla->tabla[i];
            while(actual != NULL){
                longitud++;
                actual = actual->siguiente;
            }
            
            if(longitud > 1){
                total_colisiones += (longitud - 1);
            }
            
            if(longitud > max_colisiones){
                max_colisiones = longitud;
            }
        }
    }
    
    printf("Buckets usados: %d de %d (%.2f%%)\n", 
           buckets_usados, HASH_SIZE, 
           (float)buckets_usados * 100.0 / HASH_SIZE);
    printf("Colisiones: %d\n", total_colisiones);
    printf("Max cadena: %d\n", max_colisiones);
    
    if(buckets_usados > 0){
        printf("Promedio cadena: %.2f\n", 
               (float)tabla->num_elementos / buckets_usados);
    }
    
    printf("================================\n");
}

//inicializa la tabla hash con todos los individuos del sistema
void InicializarHashTable(Mapa *grafo){
    printf("\nInicializando Hash Table...\n");
    
    if(grafo->hash_individuos != NULL){
        printf("Liberando anterior...\n");
        LiberarHashTable(grafo->hash_individuos);
    }
    
    grafo->hash_individuos = CrearHashTable();
    
    int insertados = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            Individuo *ind = territorio->individuos[i];
            if(ind != NULL){
                InsertarHash(grafo->hash_individuos, ind);
                insertados++;
            }
        }
    }
    
    printf("✓ %d individuos insertados\n", insertados);
}

void LiberarHashTable(HashTable *tabla){
    for(int i = 0; i < HASH_SIZE; i++){
        NodoHash *actual = tabla->tabla[i];
        while(actual != NULL){
            NodoHash *siguiente = actual->siguiente;
            free(actual);
            actual = siguiente;
        }
    }
    free(tabla);
}

void PruebaRendimientoHash(Mapa *grafo){
    printf("\n========== PRUEBA RENDIMIENTO ==========\n");
    
    if(grafo->hash_individuos == NULL){
        printf("Hash Table no inicializada.\n");
        return;
    }
    
    int total_individuos = 0;
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        total_individuos += grafo->territorios[t].num_individuos;
    }
    
    if(total_individuos == 0){
        printf("No hay individuos.\n");
        return;
    }
    
    int *ids_validos = (int*)malloc(total_individuos * sizeof(int));
    int idx = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL){
                ids_validos[idx++] = territorio->individuos[i]->ID;
            }
        }
    }
    
    int num_pruebas = (total_individuos < 1000) ? total_individuos : 1000;
    
    printf("Búsquedas: %d\n", num_pruebas);
    
    int encontrados = 0;
    
    for(int i = 0; i < num_pruebas; i++){
        int idx_aleatorio = (int)(Azar(0, idx - 1));
        if(idx_aleatorio >= idx) idx_aleatorio = idx - 1;
        if(idx_aleatorio < 0) idx_aleatorio = 0;
        
        int id_buscar = ids_validos[idx_aleatorio];
        
        Individuo *resultado = BuscarHash(grafo->hash_individuos, id_buscar);
        if(resultado != NULL){
            encontrados++;
        }
    }
    
    printf("\nResultados:\n");
    printf("Encontrados: %d\n", encontrados);
    printf("Tasa éxito: %.2f%%\n", (float)encontrados * 100.0 / num_pruebas);
    
    printf("\nEficiencia:\n");
    printf("Hash Table: O(1)\n");
    printf("Búsqueda lineal: O(%d)\n", total_individuos);
    printf("Mejora: ~%dx\n", total_individuos / 10);
    
    printf("\nEjemplos (primeros 5):\n");
    int mostrados = 0;
    
    for(int i = 0; i < num_pruebas && mostrados < 5; i++){
        int idx_aleatorio = (int)(Azar(0, idx - 1));
        if(idx_aleatorio >= idx) idx_aleatorio = idx - 1;
        if(idx_aleatorio < 0) idx_aleatorio = 0;
        
        int id_buscar = ids_validos[idx_aleatorio];
        Individuo *resultado = BuscarHash(grafo->hash_individuos, id_buscar);
        
        if(resultado != NULL){
            printf("  ID %d: %s (%s)\n",
                   resultado->ID,
                   resultado->Nombre,
                   grafo->territorios[resultado->Territorio_ID].Nombre);
            mostrados++;
        }
    }
    
    printf("========================================\n");
    
    free(ids_validos);
}

//=============================================================
//rutas criticas de contagio con dijkstra - o((n+m) log n)
//=============================================================

//reconstruye y muestra la ruta desde origen hasta destino
void ReconstruirRuta(Mapa *grafo, int origen, int destino, int *padre){
    if(destino == origen){
        printf("%s", grafo->territorios[origen].Nombre);
        return;
    }
    
    if(padre[destino] == -1){
        printf("(No hay ruta)");
        return;
    }
    
    ReconstruirRuta(grafo, origen, padre[destino], padre);
    printf(" -> %s", grafo->territorios[destino].Nombre);
}

void MostrarRutasCriticas(Mapa *grafo, int origen, float *distancia, int *padre){
    printf("\n========== RUTAS CRÍTICAS DESDE %s ==========\n", 
           grafo->territorios[origen].Nombre);
    
    printf("\n%-20s %10s %s\n", "Destino", "Distancia", "Ruta");
    printf("-----------------------------------------------------------\n");
    
    // Ordenar destinos por distancia
    typedef struct {
        int territorio;
        float distancia;
    } DistanciaTerritorio;
    
    DistanciaTerritorio destinos[NUM_TERRITORIOS];
    int num_destinos = 0;
    
    for(int v = 0; v < NUM_TERRITORIOS; v++){
        if(v != origen && distancia[v] < 999999.0){
            destinos[num_destinos].territorio = v;
            destinos[num_destinos].distancia = distancia[v];
            num_destinos++;
        }
    }
    
    // Ordenar por distancia (burbuja simple)
    for(int i = 0; i < num_destinos - 1; i++){
        for(int j = 0; j < num_destinos - i - 1; j++){
            if(destinos[j].distancia > destinos[j+1].distancia){
                DistanciaTerritorio temp = destinos[j];
                destinos[j] = destinos[j+1];
                destinos[j+1] = temp;
            }
        }
    }
    
    // Mostrar todos los destinos alcanzables
    for(int i = 0; i < num_destinos; i++){
        int v = destinos[i].territorio;
        printf("%-20s %10.3f  ", 
               grafo->territorios[v].Nombre, 
               distancia[v]);
        ReconstruirRuta(grafo, origen, v, padre);
        printf("\n");
    }
    
    // Identificar territorios inalcanzables
    int inalcanzables = 0;
    for(int v = 0; v < NUM_TERRITORIOS; v++){
        if(v != origen && distancia[v] >= 999999.0){
            inalcanzables++;
        }
    }
    
    if(inalcanzables > 0){
        printf("\n⚠ Territorios inalcanzables: %d\n", inalcanzables);
    }
    
    printf("\n============================================================\n");
    
    // Análisis de rutas críticas
    printf("\n--- ANÁLISIS DE RUTAS CRÍTICAS ---\n");
    
    // Ruta más corta
    if(num_destinos > 0){
        int mas_cercano = destinos[0].territorio;
        printf("\nTerritorio más cercano:\n");
        printf("  %s (distancia: %.3f)\n", 
               grafo->territorios[mas_cercano].Nombre,
               distancia[mas_cercano]);
        printf("  Ruta: ");
        ReconstruirRuta(grafo, origen, mas_cercano, padre);
        printf("\n");
    }
    
    // Ruta más larga
    if(num_destinos > 0){
        int mas_lejano = destinos[num_destinos-1].territorio;
        printf("\nTerritorio más lejano:\n");
        printf("  %s (distancia: %.3f)\n", 
               grafo->territorios[mas_lejano].Nombre,
               distancia[mas_lejano]);
        printf("  Ruta: ");
        ReconstruirRuta(grafo, origen, mas_lejano, padre);
        printf("\n");
    }
    
    // Rutas de alto riesgo (distancia corta)
    printf("\nRutas de ALTO RIESGO (distancia < 2.0):\n");
    int rutas_alto_riesgo = 0;
    for(int i = 0; i < num_destinos && destinos[i].distancia < 2.0; i++){
        int v = destinos[i].territorio;
        printf("  %d. %s (%.3f): ", 
               rutas_alto_riesgo + 1,
               grafo->territorios[v].Nombre,
               distancia[v]);
        ReconstruirRuta(grafo, origen, v, padre);
        printf("\n");
        rutas_alto_riesgo++;
        
        if(rutas_alto_riesgo >= 5) break; // Mostrar solo top 5
    }
    
    if(rutas_alto_riesgo == 0){
        printf("  (Ninguna)\n");
    }
    
    printf("\n===================================\n");
}

//encuentra las rutas mas cortas desde un territorio usando dijkstra
void AlgoritmoDijkstra(Mapa *grafo, int territorio_origen){
    printf("\n========== ALGORITMO DE DIJKSTRA ==========\n");
    printf("Calculando rutas desde: %s\n", grafo->territorios[territorio_origen].Nombre);
    
    if(territorio_origen < 0 || territorio_origen >= NUM_TERRITORIOS){
        printf("Error: Territorio inválido\n");
        return;
    }
    
    // Inicialización
    float *distancia = (float*)malloc(NUM_TERRITORIOS * sizeof(float));
    int *padre = (int*)malloc(NUM_TERRITORIOS * sizeof(int));
    int *visitado = (int*)malloc(NUM_TERRITORIOS * sizeof(int));
    
    for(int i = 0; i < NUM_TERRITORIOS; i++){
        distancia[i] = 999999.0; // Infinito
        padre[i] = -1;
        visitado[i] = 0;
    }
    
    distancia[territorio_origen] = 0.0;
    
    // Crear MinHeap
    MinHeap *heap = CrearMinHeap(NUM_TERRITORIOS);
    
    for(int v = 0; v < NUM_TERRITORIOS; v++){
        InsertarHeap(heap, v, distancia[v]);
    }
    
    // Algoritmo de Dijkstra
    int vertices_procesados = 0;
    
    while(!EstaVacio(heap) && vertices_procesados < NUM_TERRITORIOS){
        NodoHeap min_nodo = ExtraerMin(heap);
        int u = min_nodo.vertice;
        
        if(u == -1 || distancia[u] >= 999999.0) break;
        
        visitado[u] = 1;
        vertices_procesados++;
        
        // Relajar aristas adyacentes
        for(int v = 0; v < NUM_TERRITORIOS; v++){
            if(!visitado[v] && grafo->matrix[u][v] > 0.0){
                // Convertir proximidad a distancia
                // Mayor proximidad = menor distancia
                float peso_arista = 1.0 / grafo->matrix[u][v];
                
                float nueva_distancia = distancia[u] + peso_arista;
                
                if(nueva_distancia < distancia[v]){
                    distancia[v] = nueva_distancia;
                    padre[v] = u;
                    DisminuirClave(heap, v, nueva_distancia);
                }
            }
        }
    }
    
    // Mostrar resultados
    MostrarRutasCriticas(grafo, territorio_origen, distancia, padre);
    
    // Liberar memoria
    LiberarHeap(heap);
    free(distancia);
    free(padre);
    free(visitado);
}

//=============================================================
//clustering de cepas usando trie - o(n*l) construccion, o(l) busqueda
//=============================================================

int CharAIndice(char c){
    if(c >= 'a' && c <= 'z'){
        return c - 'a';
    }
    if(c >= 'A' && c <= 'Z'){
        return c - 'A';
    }
    return -1; // Caracter no válido
}

//crea un nuevo nodo vacio para el trie
NodoTrie* CrearNodoTrie(){
    NodoTrie *nodo = (NodoTrie*)malloc(sizeof(NodoTrie));
    nodo->es_final = 0;
    nodo->cepa_id = -1;
    
    for(int i = 0; i < ALPHABET_SIZE; i++){
        nodo->hijos[i] = NULL;
    }
    
    return nodo;
}

//crea un nuevo trie vacio
Trie* CrearTrie(){
    Trie *trie = (Trie*)malloc(sizeof(Trie));
    trie->raiz = CrearNodoTrie();
    trie->num_palabras = 0;
    return trie;
}

//inserta una palabra en el trie asociada a un id de cepa
void InsertarEnTrie(Trie *trie, const char *palabra, int cepa_id){
    NodoTrie *actual = trie->raiz;
    
    for(int i = 0; palabra[i] != '\0'; i++){
        int indice = CharAIndice(palabra[i]);
        
        // Saltar caracteres especiales (-, números, espacios)
        if(indice == -1) continue;
        
        if(actual->hijos[indice] == NULL){
            actual->hijos[indice] = CrearNodoTrie();
        }
        
        actual = actual->hijos[indice];
    }
    
    actual->es_final = 1;
    actual->cepa_id = cepa_id;
    trie->num_palabras++;
}

//busca una palabra en el trie y retorna el id de cepa si existe
int BuscarEnTrie(Trie *trie, const char *palabra){
    NodoTrie *actual = trie->raiz;
    
    for(int i = 0; palabra[i] != '\0'; i++){
        int indice = CharAIndice(palabra[i]);
        
        if(indice == -1) continue;
        
        if(actual->hijos[indice] == NULL){
            return -1; // No encontrado
        }
        
        actual = actual->hijos[indice];
    }
    
    if(actual != NULL && actual->es_final){
        return actual->cepa_id;
    }
    
    return -1;
}

void BuscarPorPrefijoRecursivo(NodoTrie *nodo, char *prefijo_actual, int nivel, 
                                Mapa *grafo, int *resultados, int *num_resultados){
    if(nodo == NULL) return;
    
    if(nodo->es_final && *num_resultados < 50){
        resultados[*num_resultados] = nodo->cepa_id;
        (*num_resultados)++;
    }
    
    for(int i = 0; i < ALPHABET_SIZE; i++){
        if(nodo->hijos[i] != NULL){
            prefijo_actual[nivel] = 'a' + i;
            prefijo_actual[nivel + 1] = '\0';
            
            BuscarPorPrefijoRecursivo(nodo->hijos[i], prefijo_actual, 
                                     nivel + 1, grafo, resultados, num_resultados);
        }
    }
}

//busca todas las cepas que empiezan con un prefijo dado
void BuscarPorPrefijo(Trie *trie, const char *prefijo, Mapa *grafo){
    printf("\n========== BÚSQUEDA POR PREFIJO: \"%s\" ==========\n", prefijo);
    
    NodoTrie *actual = trie->raiz;
    
    // Navegar hasta el final del prefijo
    for(int i = 0; prefijo[i] != '\0'; i++){
        int indice = CharAIndice(prefijo[i]);
        
        if(indice == -1) continue;
        
        if(actual->hijos[indice] == NULL){
            printf("\nNo se encontraron cepas con el prefijo \"%s\"\n", prefijo);
            printf("====================================================\n");
            return;
        }
        
        actual = actual->hijos[indice];
    }
    
    // Recolectar todos los resultados desde este nodo
    int resultados[50];
    int num_resultados = 0;
    char prefijo_actual[100];
    strcpy(prefijo_actual, prefijo);
    
    BuscarPorPrefijoRecursivo(actual, prefijo_actual, strlen(prefijo), 
                             grafo, resultados, &num_resultados);
    
    if(num_resultados == 0){
        printf("\nNo se encontraron cepas.\n");
    } else {
        printf("\nCepas encontradas: %d\n\n", num_resultados);
        printf("%-5s %-25s %10s %10s\n", "ID", "Nombre", "Tasa Cont.", "Tasa Mort.");
        printf("----------------------------------------------------------------\n");
        
        for(int i = 0; i < num_resultados; i++){
            int cepa_id = resultados[i];
            if(cepa_id >= 0 && cepa_id < NUM_CEPAS){
                Cepa *cepa = &grafo->cepas[cepa_id];
                printf("%-5d %-25s %10.2f %10.2f\n", 
                       cepa->ID,
                       cepa->Nombre,
                       cepa->Tasa_contagio,
                       cepa->Tasa_mortalidad);
            }
        }
    }
    
    printf("====================================================\n");
}

//inicializa el trie con todas las cepas del sistema
void InicializarTrie(Mapa *grafo){
    printf("\nInicializando Trie...\n");
    
    if(grafo->trie_cepas != NULL){
        printf("Liberando Trie anterior...\n");
        LiberarTrie(grafo->trie_cepas->raiz);
        free(grafo->trie_cepas);
    }
    
    grafo->trie_cepas = CrearTrie();
    
    // Insertar todas las cepas en el Trie
    for(int i = 0; i < NUM_CEPAS; i++){
        InsertarEnTrie(grafo->trie_cepas, grafo->cepas[i].Nombre, i);
    }
    
    printf("✓ Trie inicializado con %d cepas\n", grafo->trie_cepas->num_palabras);
}

void MostrarCepasPorCluster(Mapa *grafo){
    printf("\n========== CLUSTERING DE CEPAS ==========\n");
    
    // Agrupar cepas por prefijo
    const char *prefijos[] = {"Alpha", "Beta", "Gamma", "Delta", "Epsilon", 
                              "Zeta", "Eta", "Theta", "Iota", "Kappa"};
    int num_prefijos = 10;
    
    for(int p = 0; p < num_prefijos; p++){
        printf("\n--- Cluster: %s ---\n", prefijos[p]);
        
        int count = 0;
        for(int i = 0; i < NUM_CEPAS; i++){
            // Verificar si el nombre empieza con este prefijo
            if(strncmp(grafo->cepas[i].Nombre, prefijos[p], strlen(prefijos[p])) == 0){
                printf("  %d. %s (Cont: %.2f, Mort: %.2f)\n",
                       count + 1,
                       grafo->cepas[i].Nombre,
                       grafo->cepas[i].Tasa_contagio,
                       grafo->cepas[i].Tasa_mortalidad);
                count++;
            }
        }
        
        if(count == 0){
            printf("  (Ninguna cepa en este cluster)\n");
        } else {
            printf("  Total: %d cepas\n", count);
        }
    }
    
    printf("\n=========================================\n");
    
    // Estadísticas por cluster
    printf("\n--- ESTADÍSTICAS POR CLUSTER ---\n");
    printf("%-15s %8s %12s %12s\n", "Cluster", "Cantidad", "Avg Cont.", "Avg Mort.");
    printf("-----------------------------------------------------\n");
    
    for(int p = 0; p < num_prefijos; p++){
        int count = 0;
        float suma_cont = 0.0;
        float suma_mort = 0.0;
        
        for(int i = 0; i < NUM_CEPAS; i++){
            if(strncmp(grafo->cepas[i].Nombre, prefijos[p], strlen(prefijos[p])) == 0){
                count++;
                suma_cont += grafo->cepas[i].Tasa_contagio;
                suma_mort += grafo->cepas[i].Tasa_mortalidad;
            }
        }
        
        if(count > 0){
            printf("%-15s %8d %12.3f %12.3f\n", 
                   prefijos[p], 
                   count,
                   suma_cont / count,
                   suma_mort / count);
        }
    }
    
    printf("=====================================================\n");
}

void LiberarTrie(NodoTrie *nodo){
    if(nodo == NULL) return;
    
    for(int i = 0; i < ALPHABET_SIZE; i++){
        if(nodo->hijos[i] != NULL){
            LiberarTrie(nodo->hijos[i]);
        }
    }
    
    free(nodo);
}

//=============================================================
//quicksort - o(n log n) promedio
//=============================================================

//particiona el arreglo alrededor de un pivote
int Particionar(IndividuoOrden *arr, int inicio, int fin){
    float pivote = arr[fin].valor_orden;
    int i = inicio - 1;
    
    for(int j = inicio; j < fin; j++){
        // Orden descendente (mayor primero)
        if(arr[j].valor_orden >= pivote){
            i++;
            IndividuoOrden temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    
    IndividuoOrden temp = arr[i + 1];
    arr[i + 1] = arr[fin];
    arr[fin] = temp;
    
    return i + 1;
}

//ordena el arreglo usando quicksort recursivo
void QuickSort(IndividuoOrden *arr, int inicio, int fin){
    if(inicio < fin){
        int pi = Particionar(arr, inicio, fin);
        QuickSort(arr, inicio, pi - 1);
        QuickSort(arr, pi + 1, fin);
    }
}

//ordena individuos por riesgo descendente usando quicksort
void OrdenarPorRiesgoQuick(Mapa *grafo){
    printf("\n========== QUICKSORT POR RIESGO ==========\n");
    printf("Complejidad: O(n log n) promedio, O(n²) peor caso\n\n");
    
    int total = 0;
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        total += grafo->territorios[t].num_individuos;
    }
    
    if(total == 0){
        printf("No hay individuos.\n");
        return;
    }
    
    IndividuoOrden *lista = (IndividuoOrden*)malloc(total * sizeof(IndividuoOrden));
    int idx = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL){
                lista[idx].individuo = territorio->individuos[i];
                lista[idx].valor_orden = territorio->individuos[i]->Riesgo_inicial;
                idx++;
            }
        }
    }
    
    printf("Ordenando %d individuos con QuickSort...\n", idx);
    QuickSort(lista, 0, idx - 1);
    
    printf("\n--- Top 15 MAYOR RIESGO ---\n");
    for(int i = 0; i < 15 && i < idx; i++){
        Individuo *ind = lista[i].individuo;
        printf("%3d. %s (%s) - %.4f\n",
               i+1, ind->Nombre,
               grafo->territorios[ind->Territorio_ID].Nombre,
               lista[i].valor_orden);
    }
    
    printf("\n==========================================\n");
    free(lista);
}

//=============================================================
//quicksort para IndividuoRiesgo (usado en greedy)
//=============================================================

//particiona el arreglo de riesgo alrededor de un pivote (descendente)
int ParticionarRiesgo(IndividuoRiesgo *arr, int inicio, int fin){
    float pivote = arr[fin].riesgo_calculado;
    int i = inicio - 1;
    
    for(int j = inicio; j < fin; j++){
        //orden descendente (mayor primero)
        if(arr[j].riesgo_calculado >= pivote){
            i++;
            IndividuoRiesgo temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    
    IndividuoRiesgo temp = arr[i + 1];
    arr[i + 1] = arr[fin];
    arr[fin] = temp;
    
    return i + 1;
}

//ordena el arreglo de riesgo usando quicksort
void QuickSortRiesgo(IndividuoRiesgo *arr, int inicio, int fin){
    if(inicio < fin){
        int pi = ParticionarRiesgo(arr, inicio, fin);
        QuickSortRiesgo(arr, inicio, pi - 1);
        QuickSortRiesgo(arr, pi + 1, fin);
    }
}

//=============================================================
//quicksort por grado (descendente)
//=============================================================

//particiona el arreglo por grado (descendente)
int ParticionarGrado(IndividuoOrden *arr, int inicio, int fin){
    float pivote = arr[fin].valor_orden;
    int i = inicio - 1;
    
    for(int j = inicio; j < fin; j++){
        //orden descendente (mayor primero)
        if(arr[j].valor_orden >= pivote){
            i++;
            IndividuoOrden temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    
    IndividuoOrden temp = arr[i + 1];
    arr[i + 1] = arr[fin];
    arr[fin] = temp;
    
    return i + 1;
}

//ordena por grado usando quicksort
void QuickSortGrado(IndividuoOrden *arr, int inicio, int fin){
    if(inicio < fin){
        int pi = ParticionarGrado(arr, inicio, fin);
        QuickSortGrado(arr, inicio, pi - 1);
        QuickSortGrado(arr, pi + 1, fin);
    }
}

//=============================================================
//quicksort por territorio (ascendente)
//=============================================================

//particiona el arreglo por territorio (ascendente)
int ParticionarTerritorio(IndividuoOrden *arr, int inicio, int fin){
    float pivote = arr[fin].valor_orden;
    int i = inicio - 1;
    
    for(int j = inicio; j < fin; j++){
        //orden ascendente (menor primero)
        if(arr[j].valor_orden <= pivote){
            i++;
            IndividuoOrden temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    
    IndividuoOrden temp = arr[i + 1];
    arr[i + 1] = arr[fin];
    arr[fin] = temp;
    
    return i + 1;
}

//ordena por territorio usando quicksort
void QuickSortTerritorio(IndividuoOrden *arr, int inicio, int fin){
    if(inicio < fin){
        int pi = ParticionarTerritorio(arr, inicio, fin);
        QuickSortTerritorio(arr, inicio, pi - 1);
        QuickSortTerritorio(arr, pi + 1, fin);
    }
}

//=============================================================
//quicksort por tiempo de infeccion (no infectados al final)
//=============================================================

//particiona por tiempo de infeccion (ascendente, -1 al final)
int ParticionarTiempoInfeccion(IndividuoOrden *arr, int inicio, int fin){
    float pivote = arr[fin].valor_orden;
    int i = inicio - 1;
    
    for(int j = inicio; j < fin; j++){
        //si pivote es -1 (no infectado), todos van antes
        //si ambos son -1, mantener orden
        //si actual es -1, va despues del pivote
        //si pivote >= 0 y actual >= 0, orden ascendente
        int va_antes = 0;
        
        if(pivote < 0){
            //pivote no infectado, cualquier infectado va antes
            if(arr[j].valor_orden >= 0) va_antes = 1;
        } else {
            //pivote infectado
            if(arr[j].valor_orden >= 0 && arr[j].valor_orden <= pivote){
                va_antes = 1;
            }
        }
        
        if(va_antes){
            i++;
            IndividuoOrden temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    
    IndividuoOrden temp = arr[i + 1];
    arr[i + 1] = arr[fin];
    arr[fin] = temp;
    
    return i + 1;
}

//ordena por tiempo de infeccion usando quicksort
void QuickSortTiempoInfeccion(IndividuoOrden *arr, int inicio, int fin){
    if(inicio < fin){
        int pi = ParticionarTiempoInfeccion(arr, inicio, fin);
        QuickSortTiempoInfeccion(arr, inicio, pi - 1);
        QuickSortTiempoInfeccion(arr, pi + 1, fin);
    }
}

//=============================================================
//quicksort por nombre (alfabetico)
//=============================================================

//particiona por nombre alfabeticamente
int ParticionarNombre(IndividuoOrden *arr, int inicio, int fin){
    char *pivote = arr[fin].individuo->Nombre;
    int i = inicio - 1;
    
    for(int j = inicio; j < fin; j++){
        //orden ascendente alfabetico
        if(strcmp(arr[j].individuo->Nombre, pivote) <= 0){
            i++;
            IndividuoOrden temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    
    IndividuoOrden temp = arr[i + 1];
    arr[i + 1] = arr[fin];
    arr[fin] = temp;
    
    return i + 1;
}

//ordena por nombre usando quicksort
void QuickSortNombre(IndividuoOrden *arr, int inicio, int fin){
    if(inicio < fin){
        int pi = ParticionarNombre(arr, inicio, fin);
        QuickSortNombre(arr, inicio, pi - 1);
        QuickSortNombre(arr, pi + 1, fin);
    }
}

//=============================================================
//mergesort - o(n log n) garantizado
//=============================================================

//combina dos subarreglos ordenados en uno solo
void Merge(IndividuoOrden *arr, int inicio, int medio, int fin){
    int n1 = medio - inicio + 1;
    int n2 = fin - medio;
    
    IndividuoOrden *L = (IndividuoOrden*)malloc(n1 * sizeof(IndividuoOrden));
    IndividuoOrden *R = (IndividuoOrden*)malloc(n2 * sizeof(IndividuoOrden));
    
    for(int i = 0; i < n1; i++){
        L[i] = arr[inicio + i];
    }
    for(int j = 0; j < n2; j++){
        R[j] = arr[medio + 1 + j];
    }
    
    int i = 0, j = 0, k = inicio;
    
    while(i < n1 && j < n2){
        // Orden descendente (mayor primero)
        if(L[i].valor_orden >= R[j].valor_orden){
            arr[k] = L[i];
            i++;
        } else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
    
    while(i < n1){
        arr[k] = L[i];
        i++;
        k++;
    }
    
    while(j < n2){
        arr[k] = R[j];
        j++;
        k++;
    }
    
    free(L);
    free(R);
}

//ordena el arreglo usando mergesort recursivo
void MergeSort(IndividuoOrden *arr, int inicio, int fin){
    if(inicio < fin){
        int medio = inicio + (fin - inicio) / 2;
        MergeSort(arr, inicio, medio);
        MergeSort(arr, medio + 1, fin);
        Merge(arr, inicio, medio, fin);
    }
}

//ordena individuos por riesgo descendente usando mergesort
void OrdenarPorRiesgoMerge(Mapa *grafo){
    printf("\n========== MERGESORT POR RIESGO ==========\n");
    printf("Complejidad: O(n log n) garantizado\n");
    printf("Ventaja: Estable y predecible\n\n");
    
    int total = 0;
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        total += grafo->territorios[t].num_individuos;
    }
    
    if(total == 0){
        printf("No hay individuos.\n");
        return;
    }
    
    IndividuoOrden *lista = (IndividuoOrden*)malloc(total * sizeof(IndividuoOrden));
    int idx = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL){
                lista[idx].individuo = territorio->individuos[i];
                lista[idx].valor_orden = territorio->individuos[i]->Riesgo_inicial;
                idx++;
            }
        }
    }
    
    printf("Ordenando %d individuos con MergeSort...\n", idx);
    MergeSort(lista, 0, idx - 1);
    
    printf("\n--- Top 15 MAYOR RIESGO ---\n");
    for(int i = 0; i < 15 && i < idx; i++){
        Individuo *ind = lista[i].individuo;
        printf("%3d. %s (%s) - %.4f\n",
               i+1, ind->Nombre,
               grafo->territorios[ind->Territorio_ID].Nombre,
               lista[i].valor_orden);
    }
    
    printf("\n==========================================\n");
    free(lista);
}

//=============================================================
//heapsort - o(n log n) in-place
//=============================================================

//mantiene la propiedad heap en el subarbol con raiz en i
void HeapifyOrdenamiento(IndividuoOrden *arr, int n, int i){
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    
    // Para orden descendente, buscar el MENOR (min-heap)
    if(left < n && arr[left].valor_orden < arr[largest].valor_orden){
        largest = left;
    }
    
    if(right < n && arr[right].valor_orden < arr[largest].valor_orden){
        largest = right;
    }
    
    if(largest != i){
        IndividuoOrden temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;
        
        HeapifyOrdenamiento(arr, n, largest);
    }
}

//ordena el arreglo usando heapsort
void HeapSort(IndividuoOrden *arr, int n){
    // Construir heap
    for(int i = n / 2 - 1; i >= 0; i--){
        HeapifyOrdenamiento(arr, n, i);
    }
    
    // Extraer elementos del heap
    for(int i = n - 1; i > 0; i--){
        IndividuoOrden temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;
        
        HeapifyOrdenamiento(arr, i, 0);
    }
}

//ordena individuos por riesgo descendente usando heapsort
void OrdenarPorRiesgoHeap(Mapa *grafo){
    printf("\n========== HEAPSORT POR RIESGO ==========\n");
    printf("Complejidad: O(n log n) garantizado\n");
    printf("Ventaja: In-place, no requiere memoria extra\n\n");
    
    int total = 0;
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        total += grafo->territorios[t].num_individuos;
    }
    
    if(total == 0){
        printf("No hay individuos.\n");
        return;
    }
    
    IndividuoOrden *lista = (IndividuoOrden*)malloc(total * sizeof(IndividuoOrden));
    int idx = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL){
                lista[idx].individuo = territorio->individuos[i];
                lista[idx].valor_orden = territorio->individuos[i]->Riesgo_inicial;
                idx++;
            }
        }
    }
    
    printf("Ordenando %d individuos con HeapSort...\n", idx);
    HeapSort(lista, idx);
    
    printf("\n--- Top 15 MAYOR RIESGO ---\n");
    for(int i = 0; i < 15 && i < idx; i++){
        Individuo *ind = lista[i].individuo;
        printf("%3d. %s (%s) - %.4f\n",
               i+1, ind->Nombre,
               grafo->territorios[ind->Territorio_ID].Nombre,
               lista[i].valor_orden);
    }
    
    printf("\n==========================================\n");
    free(lista);
}

//=============================================================
//countingsort - o(n+k) lineal
//=============================================================

//ordena por grado usando counting sort - o(n+k)
void CountingSort(Mapa *grafo){
    printf("\n========== COUNTING SORT POR GRADO ==========\n");
    printf("Complejidad: O(n + k) - Lineal\n");
    printf("Ventaja: Muy rápido para rangos pequeños\n\n");
    
    int total = 0;
    int max_grado = 0;
    
    // Encontrar máximo grado
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL){
                total++;
                int grado = territorio->individuos[i]->Grado_inicial;
                if(grado > max_grado) max_grado = grado;
            }
        }
    }
    
    if(total == 0){
        printf("No hay individuos.\n");
        return;
    }
    
    printf("Rango de grados: 0 - %d\n", max_grado);
    printf("Total individuos: %d\n\n", total);
    
    // Crear array de conteo
    int *count = (int*)calloc(max_grado + 1, sizeof(int));
    
    // Contar frecuencias
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL){
                int grado = territorio->individuos[i]->Grado_inicial;
                count[grado]++;
            }
        }
    }
    
    // Mostrar distribución
    printf("--- Distribución de Grados ---\n");
    for(int g = max_grado; g >= 0; g--){
        if(count[g] > 0){
            printf("Grado %2d: %4d individuos ", g, count[g]);
            
            // Barra visual
            int barras = count[g] / 5;
            if(barras > 40) barras = 40;
            for(int b = 0; b < barras; b++){
                printf("█");
            }
            printf("\n");
        }
    }
    
    // Mostrar individuos de mayor grado
    printf("\n--- Individuos con MAYOR GRADO ---\n");
    int mostrados = 0;
    
    for(int g = max_grado; g >= 0 && mostrados < 15; g--){
        if(count[g] > 0){
            for(int t = 0; t < NUM_TERRITORIOS && mostrados < 15; t++){
                Territorio *territorio = &grafo->territorios[t];
                for(int i = 0; i < territorio->num_individuos && mostrados < 15; i++){
                    Individuo *ind = territorio->individuos[i];
                    if(ind != NULL && ind->Grado_inicial == g){
                        printf("%3d. %s (%s) - Grado: %d\n",
                               mostrados + 1,
                               ind->Nombre,
                               territorio->Nombre,
                               g);
                        mostrados++;
                    }
                }
            }
        }
    }
    
    printf("\n=========================================\n");
    free(count);
}

//=============================================================
//funciones de comparacion auxiliares
//=============================================================

//compara dos individuos por su nivel de riesgo
int CompararPorRiesgo(const void *a, const void *b){
    IndividuoOrden *ia = (IndividuoOrden*)a;
    IndividuoOrden *ib = (IndividuoOrden*)b;
    
    if(ib->valor_orden > ia->valor_orden) return 1;
    if(ib->valor_orden < ia->valor_orden) return -1;
    return 0;
}

//compara dos individuos por su grado de contactos
int CompararPorGrado(const void *a, const void *b){
    IndividuoOrden *ia = (IndividuoOrden*)a;
    IndividuoOrden *ib = (IndividuoOrden*)b;
    
    if(ib->valor_orden > ia->valor_orden) return 1;
    if(ib->valor_orden < ia->valor_orden) return -1;
    return 0;
}

//compara dos individuos por su territorio
int CompararPorTerritorio(const void *a, const void *b){
    IndividuoOrden *ia = (IndividuoOrden*)a;
    IndividuoOrden *ib = (IndividuoOrden*)b;
    
    if(ia->valor_orden > ib->valor_orden) return 1;
    if(ia->valor_orden < ib->valor_orden) return -1;
    return 0;
}

//=============================================================
//otras funciones de ordenamiento
//=============================================================

void OrdenarPorGrado(Mapa *grafo){
    printf("\n========== ORDENAMIENTO POR GRADO (QuickSort) ==========\n");
    
    int total = 0;
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        total += grafo->territorios[t].num_individuos;
    }
    
    if(total == 0) return;
    
    IndividuoOrden *lista = (IndividuoOrden*)malloc(total * sizeof(IndividuoOrden));
    int idx = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL){
                int num_contactos = 0;
                Contacto *c = territorio->individuos[i]->contactos;
                while(c != NULL){
                    num_contactos++;
                    c = c->sgt;
                }
                
                lista[idx].individuo = territorio->individuos[i];
                lista[idx].valor_orden = num_contactos;
                idx++;
            }
        }
    }
    
    //ordenar por grado usando quicksort propio
    QuickSortGrado(lista, 0, idx - 1);
    
    printf("\n--- Top 15 MÁS CONTACTOS ---\n");
    for(int i = 0; i < 15 && i < idx; i++){
        Individuo *ind = lista[i].individuo;
        printf("%3d. %s (%s) - %d contactos\n",
               i+1, ind->Nombre,
               grafo->territorios[ind->Territorio_ID].Nombre,
               (int)lista[i].valor_orden);
    }
    
    printf("====================================================\n");
    free(lista);
}

//ordena individuos por territorio de origen
void OrdenarPorTerritorio(Mapa *grafo){
    printf("\n========== ORDENAMIENTO POR TERRITORIO ==========\n");
    
    int total = 0;
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        total += grafo->territorios[t].num_individuos;
    }
    
    if(total == 0) return;
    
    IndividuoOrden *lista = (IndividuoOrden*)malloc(total * sizeof(IndividuoOrden));
    int idx = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL){
                lista[idx].individuo = territorio->individuos[i];
                lista[idx].valor_orden = territorio->individuos[i]->Territorio_ID;
                idx++;
            }
        }
    }
    
    //ordenar por territorio usando quicksort propio
    QuickSortTerritorio(lista, 0, idx - 1);
    
    printf("\nPrimeros 5 por territorio:\n");
    int territorio_actual = -1;
    int contador = 0;
    
    for(int i = 0; i < idx; i++){
        Individuo *ind = lista[i].individuo;
        
        if(ind->Territorio_ID != territorio_actual){
            territorio_actual = ind->Territorio_ID;
            contador = 0;
            printf("\n--- %s ---\n", grafo->territorios[territorio_actual].Nombre);
        }
        
        if(contador < 5){
            printf("  %s (Riesgo: %.3f)\n", ind->Nombre, ind->Riesgo_inicial);
            contador++;
        }
    }
    
    printf("\n================================================\n");
    free(lista);
}

//=============================================================
//menu de seleccion de algoritmos de ordenamiento
//=============================================================

//muestra el menu de opciones de ordenamiento
void MenuOrdenamiento(Mapa *grafo){
    printf("\n========== ALGORITMOS DE ORDENAMIENTO ==========\n");
    printf("1. QuickSort por riesgo DESC - O(n log n) promedio\n");
    printf("2. MergeSort por riesgo DESC - O(n log n) garantizado\n");
    printf("3. HeapSort por riesgo DESC - O(n log n) in-place\n");
    printf("4. Counting Sort por grado - O(n + k) lineal\n");
    printf("5. Ordenar por tiempo de infección ASC\n");
    printf("6. Ordenar por nombre ASC\n");
    printf("7. Ordenar por territorio\n");
    printf("8. Análisis estadístico\n");
    printf("9. Mostrar semillas iniciales\n");
    printf("Seleccione: ");
}

//muestra estadisticas generales del sistema
void AnalisisDatos(Mapa *grafo){
    printf("\n========== ANÁLISIS ESTADÍSTICO ==========\n");
    
    int total_individuos = 0;
    int total_infectados = 0;
    int total_recuperados = 0;
    int total_fallecidos = 0;
    float suma_riesgo = 0.0;
    int suma_grado = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            Individuo *ind = territorio->individuos[i];
            if(ind != NULL){
                total_individuos++;
                suma_riesgo += ind->Riesgo_inicial;
                suma_grado += ind->Grado_inicial;
                
                if(ind->Infectado) total_infectados++;
                if(ind->Recuperado) total_recuperados++;
                if(ind->Fallecido) total_fallecidos++;
            }
        }
    }
    
    printf("Total individuos: %d\n", total_individuos);
    printf("Infectados activos: %d (%.2f%%)\n", total_infectados, 
           total_individuos > 0 ? (float)total_infectados/total_individuos*100 : 0);
    printf("Recuperados: %d (%.2f%%)\n", total_recuperados,
           total_individuos > 0 ? (float)total_recuperados/total_individuos*100 : 0);
    printf("Fallecidos: %d (%.2f%%)\n", total_fallecidos,
           total_individuos > 0 ? (float)total_fallecidos/total_individuos*100 : 0);
    printf("Riesgo promedio: %.4f\n", total_individuos > 0 ? suma_riesgo/total_individuos : 0);
    printf("Grado promedio: %.2f\n", total_individuos > 0 ? (float)suma_grado/total_individuos : 0);
    
    printf("\n--- Por Territorio ---\n");
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        if(territorio->num_individuos > 0){
            printf("%s: %d individuos\n", territorio->Nombre, territorio->num_individuos);
        }
    }
    
    printf("\n==========================================\n");
}

//=============================================================
//funciones para las 10 semillas iniciales de contagio
//=============================================================

// FUNCION PARA INICIALIZAR LAS 10 SEMILLAS PRECARGADAS
//configura las 10 semillas iniciales de contagio precargadas
void InicializarSemillas(Mapa *grafo){
    // Distribuir semillas en diferentes territorios y cepas
    int territorios_semilla[] = {CHINA, CHINA, ITALIA, ITALIA, EUA, 
                                  ESPANA, FRANCIA, ALEMANIA, JAPON, KOREA};
    int cepas_semilla[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    
    grafo->num_semillas = 10;
    
    for(int i = 0; i < 10; i++){
        Territorio *terr = &grafo->territorios[territorios_semilla[i]];
        
        // Buscar un individuo válido en el territorio
        int individuo_idx = i % terr->num_individuos;
        if(terr->individuos[individuo_idx] != NULL){
            grafo->semillas[i].individuo_id = terr->individuos[individuo_idx]->ID;
            grafo->semillas[i].t0 = 0;  // Tiempo inicial
            grafo->semillas[i].cepa_id = cepas_semilla[i];
        }
    }
}

// FUNCION PARA APLICAR LAS SEMILLAS (infectar a los individuos)
//aplica las semillas infectando a los individuos correspondientes
void AplicarSemillas(Mapa *grafo){
    for(int i = 0; i < grafo->num_semillas; i++){
        // Buscar el individuo por ID
        for(int t = 0; t < NUM_TERRITORIOS; t++){
            Territorio *terr = &grafo->territorios[t];
            for(int j = 0; j < terr->num_individuos; j++){
                if(terr->individuos[j] != NULL && 
                   terr->individuos[j]->ID == grafo->semillas[i].individuo_id){
                    terr->individuos[j]->Infectado = 1;
                    terr->individuos[j]->t_infeccion = grafo->semillas[i].t0;
                    terr->individuos[j]->Cepa_ID = grafo->semillas[i].cepa_id;
                }
            }
        }
    }
}

// FUNCION PARA MOSTRAR LAS SEMILLAS INICIALES
//muestra la tabla de las 10 semillas iniciales
void MostrarSemillas(Mapa *grafo){
    printf("\n========== SEMILLAS INICIALES DE CONTAGIO ==========\n");
    printf("Total: %d pacientes infectados inicialmente\n\n", grafo->num_semillas);
    
    printf("%-5s %-20s %-15s %-20s %-5s\n", "No.", "Individuo", "Territorio", "Cepa", "t0");
    printf("-------------------------------------------------------------------\n");
    
    for(int i = 0; i < grafo->num_semillas; i++){
        int ind_id = grafo->semillas[i].individuo_id;
        
        // Buscar información del individuo
        for(int t = 0; t < NUM_TERRITORIOS; t++){
            Territorio *terr = &grafo->territorios[t];
            for(int j = 0; j < terr->num_individuos; j++){
                if(terr->individuos[j] != NULL && terr->individuos[j]->ID == ind_id){
                    printf("%-5d %-20s %-15s %-20s %-5d\n",
                           i + 1,
                           terr->individuos[j]->Nombre,
                           terr->Nombre,
                           grafo->cepas[grafo->semillas[i].cepa_id].Nombre,
                           grafo->semillas[i].t0);
                }
            }
        }
    }
    
    printf("===================================================\n");
}

//=============================================================
//tabla hash para busqueda de cepas en o(1)
//=============================================================

// FUNCION PARA CREAR HASH TABLE DE CEPAS
//crea una tabla hash vacia para cepas
HashTableCepas* CrearHashTableCepas(){
    HashTableCepas *tabla = (HashTableCepas*)malloc(sizeof(HashTableCepas));
    tabla->num_elementos = 0;
    
    for(int i = 0; i < NUM_CEPAS * 2; i++){
        tabla->tabla[i] = NULL;
    }
    
    return tabla;
}

// FUNCION HASH PARA CEPAS
//calcula el indice hash para un id de cepa
int FuncionHashCepa(int id){
    return id % (NUM_CEPAS * 2);
}

// FUNCION PARA INSERTAR CEPA EN HASH TABLE
//inserta una cepa en la tabla hash de cepas
void InsertarHashCepa(HashTableCepas *tabla, Cepa *cepa){
    if(tabla == NULL || cepa == NULL) return;
    
    int indice = FuncionHashCepa(cepa->ID);
    
    NodoHashCepa *nuevo = (NodoHashCepa*)malloc(sizeof(NodoHashCepa));
    nuevo->ID = cepa->ID;
    nuevo->cepa = cepa;
    nuevo->siguiente = tabla->tabla[indice];
    tabla->tabla[indice] = nuevo;
    tabla->num_elementos++;
}

// FUNCION PARA BUSCAR CEPA POR ID - O(1)
//busca una cepa por id en la tabla hash - o(1)
Cepa* BuscarHashCepa(HashTableCepas *tabla, int id){
    if(tabla == NULL) return NULL;
    
    int indice = FuncionHashCepa(id);
    
    NodoHashCepa *actual = tabla->tabla[indice];
    while(actual != NULL){
        if(actual->ID == id){
            return actual->cepa;
        }
        actual = actual->siguiente;
    }
    
    return NULL;
}

// FUNCION PARA INICIALIZAR HASH TABLE DE CEPAS
//inicializa la tabla hash con todas las cepas
void InicializarHashTableCepas(Mapa *grafo){
    printf("\nInicializando Hash Table de Cepas...\n");
    
    if(grafo->hash_cepas != NULL){
        printf("Hash Table de Cepas ya inicializada.\n");
        return;
    }
    
    grafo->hash_cepas = CrearHashTableCepas();
    
    for(int i = 0; i < grafo->num_cepas; i++){
        InsertarHashCepa(grafo->hash_cepas, &grafo->cepas[i]);
    }
    
    printf("✓ Hash Table de Cepas inicializada con %d cepas\n", grafo->hash_cepas->num_elementos);
    printf("  Complejidad de búsqueda: O(1)\n");
}

//=============================================================
//ordenamiento por tiempo de infeccion asc
//=============================================================

// FUNCION DE COMPARACION POR TIEMPO DE INFECCION
// Los no infectados (t_infeccion = -1) van al final
//compara dos individuos por tiempo de infeccion
int CompararPorTiempoInfeccion(const void *a, const void *b){
    IndividuoOrden *ia = (IndividuoOrden*)a;
    IndividuoOrden *ib = (IndividuoOrden*)b;
    
    // Si ambos no están infectados, mantener orden
    if(ia->valor_orden < 0 && ib->valor_orden < 0) return 0;
    
    // Los no infectados (-1) van al final
    if(ia->valor_orden < 0) return 1;
    if(ib->valor_orden < 0) return -1;
    
    // Ordenar ASC por tiempo de infección
    if(ia->valor_orden < ib->valor_orden) return -1;
    if(ia->valor_orden > ib->valor_orden) return 1;
    return 0;
}

// FUNCION PARA ORDENAR POR TIEMPO DE INFECCION ASC
//ordena individuos por tiempo de infeccion ascendente
void OrdenarPorTiempoInfeccion(Mapa *grafo){
    printf("\n========== ORDENAMIENTO POR TIEMPO DE INFECCIÓN ASC ==========\n");
    printf("(No infectados al final)\n\n");
    
    int total = 0;
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        total += grafo->territorios[t].num_individuos;
    }
    
    if(total == 0) return;
    
    IndividuoOrden *lista = (IndividuoOrden*)malloc(total * sizeof(IndividuoOrden));
    int idx = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL){
                lista[idx].individuo = territorio->individuos[i];
                lista[idx].valor_orden = territorio->individuos[i]->t_infeccion;
                idx++;
            }
        }
    }
    
    //ordenar por tiempo de infeccion usando quicksort propio
    QuickSortTiempoInfeccion(lista, 0, idx - 1);
    
    printf("--- Infectados ordenados por tiempo de infección ---\n");
    int mostrados = 0;
    for(int i = 0; i < idx && mostrados < 20; i++){
        Individuo *ind = lista[i].individuo;
        if(ind->t_infeccion >= 0){
            printf("%3d. %-20s | t_infección: %3d | Cepa: %s\n",
                   mostrados + 1,
                   ind->Nombre,
                   ind->t_infeccion,
                   ind->Cepa_ID >= 0 ? grafo->cepas[ind->Cepa_ID].Nombre : "N/A");
            mostrados++;
        }
    }
    
    if(mostrados == 0){
        printf("No hay individuos infectados.\n");
    }
    
    // Contar no infectados
    int no_infectados = 0;
    for(int i = 0; i < idx; i++){
        if(lista[i].individuo->t_infeccion < 0) no_infectados++;
    }
    printf("\n... %d individuos no infectados (al final)\n", no_infectados);
    
    printf("\n=============================================================\n");
    free(lista);
}

//=============================================================
//ordenamiento por nombre alfabeticamente
//=============================================================

// FUNCION DE COMPARACION POR NOMBRE
//compara dos individuos por su nombre alfabeticamente
int CompararPorNombre(const void *a, const void *b){
    IndividuoOrden *ia = (IndividuoOrden*)a;
    IndividuoOrden *ib = (IndividuoOrden*)b;
    
    return strcmp(ia->individuo->Nombre, ib->individuo->Nombre);
}

// FUNCION PARA ORDENAR POR NOMBRE ASC
//ordena individuos por nombre alfabeticamente
void OrdenarPorNombre(Mapa *grafo){
    printf("\n========== ORDENAMIENTO POR NOMBRE ASC ==========\n");
    
    int total = 0;
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        total += grafo->territorios[t].num_individuos;
    }
    
    if(total == 0) return;
    
    IndividuoOrden *lista = (IndividuoOrden*)malloc(total * sizeof(IndividuoOrden));
    int idx = 0;
    
    for(int t = 0; t < NUM_TERRITORIOS; t++){
        Territorio *territorio = &grafo->territorios[t];
        for(int i = 0; i < territorio->num_individuos; i++){
            if(territorio->individuos[i] != NULL){
                lista[idx].individuo = territorio->individuos[i];
                lista[idx].valor_orden = 0;  // No se usa para este ordenamiento
                idx++;
            }
        }
    }
    
    printf("Ordenando %d individuos por nombre...\n", idx);
    
    //ordenar por nombre usando quicksort propio
    QuickSortNombre(lista, 0, idx - 1);
    
    printf("\n--- Primeros 25 individuos (A-Z) ---\n");
    for(int i = 0; i < 25 && i < idx; i++){
        Individuo *ind = lista[i].individuo;
        printf("%3d. %-20s | %s | Riesgo: %.3f\n",
               i + 1,
               ind->Nombre,
               grafo->territorios[ind->Territorio_ID].Nombre,
               ind->Riesgo_inicial);
    }
    
    printf("\n--- Últimos 10 individuos (Z) ---\n");
    int inicio = idx > 10 ? idx - 10 : 0;
    for(int i = inicio; i < idx; i++){
        Individuo *ind = lista[i].individuo;
        printf("%3d. %-20s | %s | Riesgo: %.3f\n",
               i + 1,
               ind->Nombre,
               grafo->territorios[ind->Territorio_ID].Nombre,
               ind->Riesgo_inicial);
    }
    
    printf("\n=================================================\n");
    free(lista);
}