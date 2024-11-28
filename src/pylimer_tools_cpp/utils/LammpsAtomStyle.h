#ifndef LAMMPS_ATOM_STYLE_H
#define LAMMPS_ATOM_STYLE_H

#include <algorithm>
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

static inline std::vector<std::string> getAtomStyleStrings()
{
    return { "none",       "angle",      "atomic", "body",     "bond",
             "charge",     "dielectric", "dipole", "dpd",      "edpd",
             "electron",   "ellipsoid",  "full",   "line",     "mdpd",
             "molecular",  "peri",       "smd",    "sph",      "sphere",
             "bpm_sphere", "spin",       "tdpd",   "template", "tri",
             "wavepacket", "hybrid" };
}
static inline std::string getAtomStyleString(AtomStyle type)
{
    return getAtomStyleStrings()[type];
}

static inline AtomStyle getAtomStyleFromString(const std::string& src)
{
    std::vector<std::string> AtomStyleString = getAtomStyleStrings();
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
