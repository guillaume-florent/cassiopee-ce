  /*------------------------------------------*/
  /* Extraction des coordonnees du pt corrige */
  /*------------------------------------------*/
  FldArrayF* coordxPC; FldArrayF* coordyPC; FldArrayF* coordzPC;
  E_Int okc1 = K_NUMPY::getFromNumpyArray(pyArrayXPC, coordxPC, true);
  E_Int okc2 = K_NUMPY::getFromNumpyArray(pyArrayYPC, coordyPC, true);
  E_Int okc3 = K_NUMPY::getFromNumpyArray(pyArrayZPC, coordzPC, true);
  E_Float* xPC = coordxPC->begin();
  E_Float* yPC = coordyPC->begin();
  E_Float* zPC = coordzPC->begin();
  
  /*----------------------------------------*/
  /* Extraction des coordonnees du pt paroi */
  /*----------------------------------------*/
  FldArrayF* coordxPW;
  FldArrayF* coordyPW;
  FldArrayF* coordzPW;
  E_Int okw1 = K_NUMPY::getFromNumpyArray(pyArrayXPW, coordxPW, true);
  E_Int okw2 = K_NUMPY::getFromNumpyArray(pyArrayYPW, coordyPW, true);
  E_Int okw3 = K_NUMPY::getFromNumpyArray(pyArrayZPW, coordzPW, true);
  E_Float* xPW = coordxPW->begin();
  E_Float* yPW = coordyPW->begin();
  E_Float* zPW = coordzPW->begin();

  /*--------------------------------------------*/
  /* Extraction des coordonnees du pt interpole */
  /*--------------------------------------------*/
  FldArrayF* coordxPI; FldArrayF* coordyPI; FldArrayF* coordzPI;
  E_Int oki1 = K_NUMPY::getFromNumpyArray(pyArrayXPI, coordxPI, true);
  E_Int oki2 = K_NUMPY::getFromNumpyArray(pyArrayYPI, coordyPI, true);
  E_Int oki3 = K_NUMPY::getFromNumpyArray(pyArrayZPI, coordzPI, true);
  E_Float* xPI = coordxPI->begin();
  E_Float* yPI = coordyPI->begin();
  E_Float* zPI = coordzPI->begin();

  /*--------------------------------------------*/
  /* Extraction des var posttraitememnt */
  /*--------------------------------------------*/
  FldArrayF* densF; FldArrayF* pressF; FldArrayF* utauF; FldArrayF* yplusF;
  E_Int okD = K_NUMPY::getFromNumpyArray(pyArrayDens    , densF , true);
  E_Int okP = K_NUMPY::getFromNumpyArray(pyArrayPressure, pressF, true);
  E_Float* density  = densF->begin();
  E_Float* pressure = pressF->begin();
  E_Int okU, okY;  E_Float* utau; E_Float* yplus;

  okU = K_NUMPY::getFromNumpyArray(pyArrayUtau , utauF , true);
  okY = K_NUMPY::getFromNumpyArray(pyArrayYplus, yplusF, true);
  utau  = utauF->begin();
  yplus = yplusF->begin();
/*  if(bcType != 0)
  {
    okU = K_NUMPY::getFromNumpyArray(pyArrayUtau , utauF , true);
    okY = K_NUMPY::getFromNumpyArray(pyArrayYplus, yplusF, true);
    utau  = utauF->begin();
    yplus = yplusF->begin();
  }
  else
  {
    utau  = pressF->begin();
    yplus = pressF->begin();
  }
*/

    
