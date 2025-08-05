#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace pylimer_tools::utils {
#define LAMMPS_ATOM_STYLES                                                     \
  X(NONE, "none")                                                              \
  X(ANGLE, "angle")                                                            \
  X(ATOMIC, "atomic")                                                          \
  X(BODY, "body")                                                              \
  X(BOND, "bond")                                                              \
  X(BPM_SPHERE, "bpm/sphere")                                                  \
  X(CHARGE, "charge")                                                          \
  X(DIELECTRIC, "dielectric")                                                  \
  X(DIPOLE, "dipole")                                                          \
  X(DPD, "dpd")                                                                \
  X(EDPD, "edpd")                                                              \
  X(ELECTRON, "electron")                                                      \
  X(ELLIPSOID, "ellipsoid")                                                    \
  X(FULL, "full")                                                              \
  X(LINE, "line")                                                              \
  X(MDPD, "mdpd")                                                              \
  X(MOLECULAR, "molecular")                                                    \
  X(PERI, "peri")                                                              \
  X(RHEO, "rheo")                                                              \
  X(RHEO_THERMAL, "rheo/thermal")                                              \
  X(SMD, "smd")                                                                \
  X(SPH, "sph")                                                                \
  X(SPHERE, "sphere")                                                          \
  X(SPIN, "spin")                                                              \
  X(TDPD, "tdpd")                                                              \
  X(TEMPLATE, "template")                                                      \
  X(TRI, "tri")                                                                \
  X(WAVEPACKET, "wavepacket")                                                  \
  X(HYBRID, "hybrid")

enum AtomStyle : int
{
#define X(e, s) e,
  LAMMPS_ATOM_STYLES
#undef X
};

static inline std::string
getAtomStyleString(const AtomStyle type)
{
  switch (type) {
#define X(e, s)                                                                \
  case e:                                                                      \
    return s;
    LAMMPS_ATOM_STYLES
#undef X
    default:
      throw std::runtime_error("Atom style " + std::to_string(type) +
                               " not recognized");
  }
}

static inline std::vector<std::string>
getAtomStyleStrings()
{
  return {
#define X(e, s) s,
    LAMMPS_ATOM_STYLES
#undef X
  };
}

static inline AtomStyle
getAtomStyleFromString(const std::string& src)
{
  // unfortunately, there's no support for switch(std::string) in C++, currently
#define X(e, s)                                                                \
  if (src == s) {                                                              \
    return e;                                                                  \
  }
  LAMMPS_ATOM_STYLES
#undef X
  throw std::invalid_argument("Atom style " + src + " not recognized");
}
}
