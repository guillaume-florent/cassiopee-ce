/*
    Copyright 2013-2017 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

//=============================================================================
// Keep information on the transcient software state
//=============================================================================
#ifndef _CPLOT_STATE_H_
#define _CPLOT_STATE_H_
#include <pthread.h>
#include <sys/time.h>
#include <cstdio>
#include "Def/DefTypes.h"

#if defined( _WIN32 ) || defined( _WIN64 )
#include <winsock.h>
#endif

/* Define a chained integer list */
struct chain_int {
    int        value;
    chain_int* next;
};

/* Define informations on the state of software */
struct CPlotState {
    int mode;  // mode=MESH, SOLID, RENDER, SCALARFIELD, VECTORFIELD
    int dim;   // 1D, 2D, 3D display

    int scalarField;                               // variable number for scalar field
    int vectorField1, vectorField2, vectorField3;  // variables for vector field

    int ijk1D;  // i,j,k mode (1D): i(0), j(1), k(2)
    int var1D;  // 1D secondary variable (inutilise)
    int var2D;  // 2D secondary variable (toggle entre xy, yz, xz, pour dim=2)

    int ghostifyDeactivatedZones;  // display deactivated zones as blended
    int edgifyDeactivatedZones;    // display deactivated zones as edges (NI)
    int edgifyActivatedZones;      // display activated zones as edges (NI)

    float alpha;            // parameter for alpha blending (1: no alpha blending)
    int   currentMenuItem;  // current displayed menu item
    int   inmenu;           // we are in menu mode or not?
    int   smoke;            // if particles smoke is to be displayed

    int render;             // -1: dont push data, dont flush buffer
                            // 0: push but dont flush
                            // 1, push and flush at end of display
    int   bb;               // if 1, display BB in global display
    int   header;           // if 1, display header in global display
    int   info;             // if 1, display the info on selected zone
    int   texture;          // if 1, display textures
    int   menu;             // if 1, display bottom menu
    int   isoLegend;        // if 1, display iso legend
    int   offscreen;        // 1 offscreen rendering with mesa, 2 with FBO
    void* offscreenBuffer;  // buffer for rendering (offscreen)

    // overlay message
    char* message;

    // Stereo
    int    stereo;      // stereo mode (anaglyph)
    double stereoDist;  // stereo distance

    // Shadow mapping
    int shadow;

    // Light position (for light and shadow)
    double lightOffsetX;
    double lightOffsetY;

    // Depth of field
    int    DOF;
    double dofPower;

    // Last selected zone
    int selectedZone;
    // List of deactivated zones
    chain_int* deactivatedZones;
    // Clicked point
    double  activePointX, activePointY, activePointZ, activePointZBuf;
    int     activePointI, activePointJ, activePointK;
    double* activePointF;
    // Pressed keys
    int  kkeysActivated;  // if 1, short cut keys are activated, 0: deactivated
    int  kcursor;
    char keys[128];  // store pressed keys

    // Textures
    char shaderPath[1028];  // chemin pour les fichiers + shaders
    char envmapFile[120];
    int  updateEnvmap;      // dit a display d'updater l'envmap
    int  updateBackground;  // dit a displau d'updater la texture de background

    // Billboard (material=9)
    int billBoardType;  // type du billboard en memoire (texture)
    int billBoardNi;    // nbre de samples en i
    int billBoardNj;    // nbre de samples en j
    int billBoardD;     // no du champ a utiliser pour la taille des billboards
    int billBoardT;     // no du champ a utiliser pour le choix du billboard

    // Autres
    char winTitle[120];        // titre de la fenetre
    char file[120];            // fichier en cours de display (si il y en a)
    char localPathName[1028];  // local directory name

    // cursor
    int cursorType;    // current cursor type
    int updateCursor;  // dit a display d'updater le curseur de la souris

    // Position de la souris
    int   activeMouseButton;  // bouton de la souris presse (dernier click)
    int   activeMouseX;       // position de la souris (au dernier click)
    int   activeMouseY;
    int   currentMouseButton;  // bouton courant
    float currentMousePosX;    // position courante de la souris
    float currentMousePosY;    // quand on est en drag mode (shift + click)
    float currentMousePosZ;
    int   modifier;  // modifieur quand la souris est pressee

    // Styles
    int   meshStyle;       // style pour le mesh mode
    int   solidStyle;      // style pour le solid mode
    int   scalarStyle;     // style pour le scalar mode
    int   vectorStyle;     // style pour le vector mode
    int   selectionStyle;  // style pour la selection (0: bleue, 1: alpha)
    int   colormap;        // colormap
    int   isoLight;        // light ou pas light -> isoLight
    int   niso;            // nbre d'isos (global)
    float isoEdges;        // edge entre les isos

    virtual ~CPlotState( ) {}

    // lock=1 pendant le display, les donnees ne doivent alors
    // pas etre modifiees!
    volatile int    lock;
    pthread_mutex_t lock_mutex;
    pthread_cond_t  unlocked_display;
    void            lockDisplay( ) {
        pthread_mutex_lock( &lock_mutex );
        lock = 1;
        pthread_mutex_unlock( &lock_mutex );
    }
    void unlockDisplay( ) {
        pthread_mutex_lock( &lock_mutex );
        lock = 0;
        pthread_cond_signal( &unlocked_display );
        pthread_mutex_unlock( &lock_mutex );
    }
    void syncDisplay( int time = 5 ) {
        pthread_mutex_lock( &lock_mutex );
        if ( lock == 1 ) pthread_cond_wait( &unlocked_display, &lock_mutex );
        pthread_mutex_unlock( &lock_mutex );
    }

    int    freeGPURes;      // dit a display de liberer les ressources du GPU
    E_LONG freeGPUData[4];  // Etats de liberation des ressources GPUs, independant de la methode de rendering
    // 0: mode, 1: start ou size (ptr), 2: end ou -1 (ptr),
    // 3: permanent
    int*            freeGPUPtr;  // ptr
    virtual void    freeGPUResources( ) = 0;
    pthread_mutex_t gpures_mutex;
    pthread_cond_t  unlocked_gpures;
    volatile int    _lockGPURes;
    // lockGPURes=display1: On bloque les ressources GPUs qui ne doivent pas etre
    //                      modifiees pendant ce temps
    void lockGPURes( ) {
        pthread_mutex_lock( &gpures_mutex );
        _lockGPURes = 1;
        pthread_mutex_unlock( &gpures_mutex );
    }
    void unlockGPURes( ) {
        pthread_mutex_lock( &gpures_mutex );
        _lockGPURes = 0;
        pthread_cond_signal( &unlocked_gpures );
        pthread_mutex_unlock( &gpures_mutex );
    }
    virtual void syncGPURes( int dt = 5 ) {
        pthread_mutex_lock( &gpures_mutex );
        if ( _lockGPURes == 1 ) pthread_cond_wait( &unlocked_gpures, &gpures_mutex );
        pthread_mutex_unlock( &gpures_mutex );
    }

    int timeStep;           // temps en ms entre 2 appels de display
    int mouseFirstClicked;  // true if mouse has been clicked already before
                            // redisplay

    // Export
    char            exportFile[120];            // nom du fichier pour l'export
    int             exportWidth, exportHeight;  // resolution de l'export
    int             continuousExport;           // =0, pas d'export continue, =1, export continue
    FILE*           ptrFile;                    // ptr fichier pour les exports continus
    void*           context;                    // contexte pour les exports continus
    int             exportNumber;               // no ajoute au nom du fichier
    int             shootScreen;                // dit a display de shooter
    pthread_mutex_t export_mutex;
    pthread_cond_t  unlocked_export;

    // Others
    int fullScreen;  // 1: full screen, 0: window
    int bgColor;     // background color
    int autoblank;   // if 1, cplot try to set automatically the
                     // blanking plugin

    // 1D overlay display size
    int gridSizeI;  // nombre de slots 1D affichables sur le display
    int gridSizeJ;

    // animated shaders
    unsigned int ktime;   // the current time, incremented at each display
    unsigned int ktimer;  // if ktimer!=0, redisplay is triggered, timer
                          // is decremented at each display
};

// Init software state
void initState( );
#endif
