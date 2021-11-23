#ifndef _MODELDEF_H
#define _MODELDEF_H

class ModelDef {
  public:
    /* Constructors */
    ModelDef() :
      flag_full_pic (true),
      flag_isothermal (true),
      gamma (1.),
      phi_ref (0.),
      n_ref (1.),
      temp_ref(1.) {} 

    /* Public methods */
    bool is_full_pic () const { return flag_full_pic; }
    bool is_isothermal () const { return !is_full_pic() && flag_isothermal; }

    /* Member */
    bool flag_full_pic;
    bool flag_isothermal;
    Real gamma;
    Real phi_ref;
    Real n_ref;
    Real temp_ref;
};

#endif
