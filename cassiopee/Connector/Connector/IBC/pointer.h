E_Float* utau_vec    = utauPtr  + ideb;
E_Float* ro_vec      = densPtr  + ideb;
E_Float* yplus_vec   = yplusPtr + ideb;
E_Float* press_vec   = pressPtr + ideb;

E_Float* aa_vec       = tmp + (ithread-1)*size*13;
E_Float* nutcible_vec = aa_vec + size;
E_Float* utauv_vec    = aa_vec + size*2;
E_Float* ucible_vec   = aa_vec + size*3;
E_Float* vcible_vec   = aa_vec + size*4;
E_Float* wcible_vec   = aa_vec + size*5;
E_Float* ut_vec       = aa_vec + size*6;
E_Float* vt_vec       = aa_vec + size*7;
E_Float* wt_vec       = aa_vec + size*8;
E_Float* sign_vec     = aa_vec + size*9;
E_Float*    mu_vec    = aa_vec + size*10;
E_Float*  uext_vec    = aa_vec + size*11;
E_Float* alpha_vec    = aa_vec + size*12;
