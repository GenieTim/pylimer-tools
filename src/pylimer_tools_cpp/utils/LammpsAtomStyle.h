#ifndef LAMMPS_ATOM_STYLE
#define LAMMPS_ATOM_STYLE

#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {
  enum AtomStyle : int
  {
    NONE,
    ANGLE,
    ATOMIC,
    BODY,
    BOND,
    CHARGE,
    DIELECTRIC,
    DIPOLE,
    DPD,
    EDPD,
    ELECTRON,
    ELLIPSOID,
    FULL,
    LINE,
    MDPD,
    MOLECULAR,
    PERI,
    SMD,
    SPH,
    SPHERE,
    BPM_SPHERE,
    SPIN,
    TDPD,
    TEMPLATE,
    TRI,
    WAVEPACKET,
    HYBRID
  };

  std::vector<std::string> AtomStyleString = {
    "none",       "angle",      "atomic", "body",     "bond",
    "charge",     "dielectric", "dipole", "dpd",      "edpd",
    "electron",   "ellipsoid",  "full",   "line",     "mdpd",
    "molecular",  "peri",       "smd",    "sph",      "sphere",
    "bpm_sphere", "spin",       "tdpd",   "template", "tri",
    "wavepacket", "hybrid"
  };

  std::string getAtomStyleString(AtomStyle type)
  {
    return AtomStyleString[type];
  }

  AtomStyle getAtomStyleFromString(const std::string &src) {
    auto it = std::find(AtomStyleString.begin(), AtomStyleString.end(), src);
    if (it != AtomStyleString.end()) {
      int index = it - AtomStyleString.begin();
      return static_cast<AtomStyle>(index);
    }
    return AtomStyle::NONE;
  }
}
}

#endif
