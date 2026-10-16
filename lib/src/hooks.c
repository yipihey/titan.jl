/* PORT phase 2: swappable routine hooks.
 * In libtitan.so the Fortran eos_/opac_ objects are renamed (objcopy) to
 * eos_f_/opac_f_; these shims dispatch to a Julia @cfunction if one is set,
 * otherwise to the original Fortran.  titan.x is linked without this file.  */
typedef void (*titan_fp)(void);
extern void eos_f_(void);
extern void opac_f_(void);
static titan_fp hook_eos = 0, hook_opac = 0;
void titan_set_hook(int which, titan_fp f) {
    if (which == 1) hook_eos = f;
    else if (which == 2) hook_opac = f;
}
void eos_(void)  { if (hook_eos)  hook_eos();  else eos_f_();  }
void opac_(void) { if (hook_opac) hook_opac(); else opac_f_(); }
