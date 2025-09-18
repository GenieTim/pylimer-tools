#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/io/AveFileReader.h"
#include "../../src/pylimer_tools_cpp/io/DataFileParser.h"
#include "../../src/pylimer_tools_cpp/io/DataFileWriter.h"
#include "../../src/pylimer_tools_cpp/io/DumpFileParser.h"
#include "catch2/matchers/catch_matchers_container_properties.hpp"
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>
extern "C"
{
#include <igraph/igraph.h>
}

namespace pu = pylimer_tools::utils;
namespace pe = pylimer_tools::entities;

TEST_CASE("FileParsers can be used", "[utils][DumpFileParser][DataFileParser]")
{
  std::cout << "Running test \"FileParsers can be used\"" << std::endl;
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  CHECK(std::filesystem::exists(suspectedPath));

  SECTION("Reading from dump file works")
  {
    std::string inputFile = suspectedPath + "/lammps_dump_small.lammpstrj";
    pu::DumpFileParser parser = pu::DumpFileParser(inputFile);
    CHECK(parser.getFilePath() == inputFile);
    CHECK_THROWS(parser.hasKey("BOX BOUNDS"));
    CHECK(parser.getLength() == 1);
    CHECK_NOTHROW(parser.read());
    CHECK(parser.getLength() == 1);
    CHECK_THROWS(parser.readNGroups(1, 10));
    CHECK_THROWS(parser.readNGroups(10, 1));
    CHECK_THROWS(parser.getValuesForAt<int>(0, "BOX BOUNDS", "notexisting"));
    // call copy constructor
    pu::DumpFileParser parser2 = parser;
    CHECK(parser2.getLength() == 1);

    // test other reading capabilities
    pu::DumpFileParser parser3 =
      pu::DumpFileParser(suspectedPath + "/lammps_dump_small.lammpstrj");
    CHECK(parser3.getValuesForAt<double>(0, "BOX BOUNDS", 1).size() == 3);
    CHECK_THROWS(parser3.getValuesForAt<double>(0, "NOT EXISTING", 9));
    std::vector<long int> timeSteps = parser3.readTimeSteps();
    CHECK(timeSteps.size() == 1);
    CHECK(timeSteps[0] == 70764);
    std::vector<std::vector<pe::Atom>> atoms = parser3.readAtoms();
    CHECK(atoms.size() == 1);
    CHECK(atoms[0].size() == 12);
    CHECK(atoms[0][2].getId() == 30000);
    std::vector<pe::Box> boxes = parser3.readBoxes();
    CHECK(boxes.size() == 1);
    CHECK_THAT(boxes[0].getLx(),
               Catch::Matchers::WithinRel(4.8545999999999999e+01));

    // test throws
    CHECK_THROWS(pu::DumpFileParser("not-existing-file.out"));

    // test without atoms
    pu::DumpFileParser parser4 = pu::DumpFileParser(
      suspectedPath + "/lammps_dump_small_no_atoms.lammpstrj");
    CHECK(parser4.getLength() == 1);
    CHECK(parser.getValuesForAt<double>(0, "BOX BOUNDS", 1).size() == 3);
    CHECK(parser.getValuesForAt<double>(0, "BOX BOUNDS", 1)[0] ==
          4.8545999999999999e+01);
  }

  SECTION("Reading from data files works")
  {
    pu::DataFileParser parser = pu::DataFileParser();
    parser.read(suspectedPath + "/lammps_data_file.out");
    CHECK(parser.getNrOfAtoms() == 3000);

    // call copy constructor
    pu::DataFileParser parser2 = pu::DataFileParser();
    parser2 = parser;
    CHECK(parser2.getNrOfAtoms() == 3000);
    // test throws
    CHECK_THROWS(parser.read("not-existing-file.out"));

    // test angles reading
    pu::DataFileParser parser3 = pu::DataFileParser();
    parser3.read(suspectedPath + "/lammps_data_file_small_wangles.out");
    CHECK(parser3.getNrOfAngles() == 1);

    // BENCHMARK("DataFileParserOld")
    // {
    //   pu::DataFileParser parser4 = pu::DataFileParser();
    //   parser4.read(suspectedPath + "/big_dump_file_data.out");
    //   return parser4.getNrOfAngles();
    // };

    // BENCHMARK("DataFileParserNew")
    // {
    //   pu::DataFileParser2 parser5 = pu::DataFileParser2();
    //   parser5.read(suspectedPath + "/big_dump_file_data.out");
    //   return parser5.getNrOfAngles();
    // };
  }

  SECTION("Reading large files is sensibly fast")
  {
    pu::DumpFileParser parser =
      pu::DumpFileParser(suspectedPath + "/big_dump_file.lammpstrj");
    // pre-read multiple
    CHECK_NOTHROW(parser.readNGroups(9, 12));
    CHECK(parser.hasKey("BOX BOUNDS") == true);
    CHECK(parser.hasKey("NO EXISTING") == false);

    // pu::DumpFileParser2 parser5 =
    //   pu::DumpFileParser2(suspectedPath + "/big_dump_file.lammpstrj");
    // CHECK(parser5.getLength() == parser.getLength());

    // BENCHMARK("DumpFileParserOld")
    // {
    //   pu::DumpFileParser parser4 =
    //     pu::DumpFileParser(suspectedPath + "/big_dump_file.lammpstrj");
    //   return parser4.getLength();
    // };

    // BENCHMARK("DumpFileParserNew")
    // {
    //   pu::DumpFileParser2 parser5 =
    //     pu::DumpFileParser2(suspectedPath + "/big_dump_file.lammpstrj");
    //   return parser5.getLength();
    // };
  }
};

TEST_CASE("Writers can be used", "[utils][DataFileWriter][DataFileParser]")
{
  std::cout << "Running test \"Writers can be used\"" << std::endl;
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  CHECK(std::filesystem::exists(suspectedPath));

  SECTION("Files are read and written")
  {
    // TODO: implement
    pe::UniverseSequence universeSeq = pe::UniverseSequence();
    std::string largeInputFile =
      suspectedPath + "/structure/network_100_a_46.structure.out";
    universeSeq.initializeFromDataSequence({ { largeInputFile } });
    pe::Universe universe = universeSeq.atIndex(0);

    // add angles
    auto angles = universe.detectAngles();
    std::vector<int> angleTypes;
    angleTypes.reserve(angles["angle_from"].size());
    for (size_t i = 0; i < angles["angle_from"].size(); i++) {
      angleTypes.push_back(i % 4);
    }
    universe.addAngles(angles["angle_from"],
                       angles["angle_via"],
                       angles["angle_to"],
                       angleTypes);
    CHECK(universe.getNrOfAngles() > 0);

    // add dihedrals
    auto detectedAngles = universe.detectDihedralAngles();
    CHECK(detectedAngles["dihedral_angle_from"].size() > 0);
    std::vector<int> dihedralAngleTypes;
    dihedralAngleTypes.reserve(detectedAngles["dihedral_angle_from"].size());
    for (size_t i = 0; i < detectedAngles["dihedral_angle_from"].size(); i++) {
      dihedralAngleTypes.push_back(1);
    }
    universe.addDihedralAngles(detectedAngles["dihedral_angle_from"],
                               detectedAngles["dihedral_angle_via1"],
                               detectedAngles["dihedral_angle_via2"],
                               detectedAngles["dihedral_angle_to"],
                               dihedralAngleTypes);

    universe.setPropertyValue(1, "charge", 1.05);

    // write data file
    pu::DataFileWriter writer = pu::DataFileWriter(universe);
    writer.configIncludeAngles(true);
    writer.configIncludeDihedralAngles(true);
    writer.configMoleculeIdxForSwap(false);
    writer.configMoveIntoBox(true);
    writer.configAttemptImageReset(true);
    std::string fileToWrite = suspectedPath + "/tmp_data_file.structure.out";
    writer.setCustomAtomFormat(
      "$atomId\t$atomType\t$charge\t$x\t$y\t$z\t$nx\t$ny\t$nz");
    writer.writeToFile(fileToWrite);

    pe::UniverseSequence seq = pe::UniverseSequence();
    seq.setDataFileAtomStyle({ pu::AtomStyle::CHARGE });
    seq.initializeFromDataSequence({ { fileToWrite } });
    pe::Universe readUniverse = seq.atIndex(0);

    CHECK(universe.getNrOfAtoms() == readUniverse.getNrOfAtoms());
    CHECK(universe.getNrOfBonds() == readUniverse.getNrOfBonds());
    CHECK(universe.getNrOfAngles() == readUniverse.getNrOfAngles());
    CHECK(universe.getNrOfDihedralAngles() ==
          readUniverse.getNrOfDihedralAngles());

    auto newBonds = readUniverse.getBonds();
    auto oldBonds = universe.getBonds();

    for (size_t i = 0; i < newBonds.size(); i++) {
      CHECK(newBonds["bond_type"][i] == oldBonds["bond_type"][i]);
      CHECK(newBonds["bond_from"][i] == oldBonds["bond_from"][i]);
      CHECK(newBonds["bond_to"][i] == oldBonds["bond_to"][i]);
    }

    auto newAngles = readUniverse.getAngles();
    auto oldAngles = universe.getAngles();

    for (size_t i = 0; i < newAngles.size(); i++) {
      CHECK(newAngles["angle_type"][i] == oldAngles["angle_type"][i]);
      CHECK(newAngles["angle_from"][i] == oldAngles["angle_from"][i]);
      CHECK(newAngles["angle_via"][i] == oldAngles["angle_via"][i]);
      CHECK(newAngles["angle_to"][i] == oldAngles["angle_to"][i]);
    }

    // std::filesystem::remove(fileToWrite);
  }
}

TEST_CASE("All atom types are supported",
          "[utils][DataFileWriter][DataFileParser][io]")
{
  std::cout << "Running test \"All atom types are supported\"" << std::endl;
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  CHECK(std::filesystem::exists(suspectedPath));

  std::unordered_map<pylimer_tools::utils::AtomStyle, std::vector<std::string>>
    atomPropertiesPerType = {
      { pu::AtomStyle::ANGLE, {} },  // molecule-ID already handled separately
      { pu::AtomStyle::ATOMIC, {} }, // no additional properties
      { pu::AtomStyle::BODY, { "bodyflag", "mass" } },
      { pu::AtomStyle::BOND, {} }, // molecule-ID already handled separately
      { pu::AtomStyle::BPM_SPHERE, { "diameter", "density" } },
      { pu::AtomStyle::CHARGE, { "charge" } },
      { pu::AtomStyle::DIELECTRIC,
        { "charge",
          "mux",
          "muy",
          "muz",
          "area",
          "ed",
          "em",
          "epsilon",
          "curvature" } },
      { pu::AtomStyle::DIPOLE, { "charge", "mux", "muy", "muz" } },
      { pu::AtomStyle::DPD, { "theta" } },
      { pu::AtomStyle::EDPD, { "edpd_temp", "edpd_cv" } },
      { pu::AtomStyle::ELECTRON, { "charge", "espin", "eradius" } },
      { pu::AtomStyle::ELLIPSOID, { "ellipsoidflag", "density" } },
      { pu::AtomStyle::FULL,
        { "charge" } }, // molecule-ID already handled separately
      { pu::AtomStyle::LINE, { "lineflag", "density" } },
      { pu::AtomStyle::MDPD, { "rho" } },
      { pu::AtomStyle::MOLECULAR,
        {} }, // molecule-ID already handled separately
      { pu::AtomStyle::PERI, { "volume", "density" } },
      { pu::AtomStyle::RHEO, { "status", "rho" } },
      { pu::AtomStyle::RHEO_THERMAL, { "status", "rho", "energy" } },
      { pu::AtomStyle::SMD,
        { "volume", "mass", "kradius", "cradius", "x0", "y0", "z0" } },
      { pu::AtomStyle::SPH, { "rho", "esph", "cv" } },
      { pu::AtomStyle::SPHERE, { "diameter", "density" } },
      { pu::AtomStyle::SPIN, { "spx", "spy", "spz", "sp" } },
      // TDPD has variable number of arguments, cannot test like this
      { pu::AtomStyle::TEMPLATE, { "template_index", "template_atom" } },
      { pu::AtomStyle::TRI, { "triangleflag", "density" } },
      { pu::AtomStyle::WAVEPACKET,
        { "charge", "espin", "eradius", "etag", "cs_re", "cs_im" } },
      // HYBRID is complex and typically not tested in isolation
    };

  // Test each atom style
  for (const auto& [atomStyle, properties] : atomPropertiesPerType) {
    SECTION("Testing atom style: " + pu::getAtomStyleString(atomStyle))
    {
      pe::Universe universe = pe::Universe(10., 10., 10.);

      // Add basic atoms
      universe.addAtoms({ 1, 2, 3, 4, 5 },
                        { 0, 1, 0, 2, 0 },
                        { 1., 2., 3., 4., 5. },
                        { 0., 0., 0., 0., 0. },
                        { 0., 1., 2., 3., 4. },
                        { 0, 0, 0, 1, 1 },
                        { 0, 0, 0, 0, 0 },
                        { 0, 0, 0, 0, 0 });

      // Add bonds for testing
      universe.addBonds({ 1, 2, 3, 4 }, { 2, 3, 4, 5 }, { 1, 1, 1, 1 });

      // Add required properties for this atom style
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<double> realDist(0.1, 10.0);
      std::uniform_int_distribution<int> intDist(0, 5);

      for (const std::string& property : properties) {
        for (int atomId = 1; atomId <= 5; ++atomId) {
          if (property.find("flag") != std::string::npos ||
              property == "status" || property == "template_index" ||
              property == "template_atom" || property == "etag") {
            // Integer properties
            int intVal = intDist(gen);
            universe.setPropertyValue(atomId - 1, property.c_str(), intVal);
            CHECK_THAT(universe.getAtom(atomId).getExtraData()[property],
                       Catch::Matchers::WithinAbs(intVal, 1e-9));
          } else {
            // Double properties
            int doubleVal = realDist(gen);
            universe.setPropertyValue(atomId - 1, property.c_str(), doubleVal);
            CHECK_THAT(universe.getAtom(atomId).getExtraData()[property],
                       Catch::Matchers::WithinRel(doubleVal, 1e-9));
          }
        }
      }

      // Special handling for molecule-based atom styles
      if (atomStyle == pu::AtomStyle::ANGLE ||
          atomStyle == pu::AtomStyle::BOND ||
          atomStyle == pu::AtomStyle::BPM_SPHERE ||
          atomStyle == pu::AtomStyle::FULL ||
          atomStyle == pu::AtomStyle::LINE ||
          atomStyle == pu::AtomStyle::MOLECULAR ||
          atomStyle == pu::AtomStyle::TEMPLATE ||
          atomStyle == pu::AtomStyle::TRI) {
        // These styles include molecule-ID which is handled separately
        for (int atomId = 1; atomId <= 5; ++atomId) {
          universe.setPropertyValue(atomId - 1, "molecule_id", atomId % 3);
        }
      }

      // Write data file
      pu::DataFileWriter writer = pu::DataFileWriter(universe);
      writer.configAtomStyle(atomStyle);

      std::string atomStyleString = pu::getAtomStyleString(atomStyle);
      std::ranges::replace(atomStyleString, '/', '_');
      std::string fileToWrite =
        suspectedPath + "/tmp_data_file_" + atomStyleString + ".structure.out";

      writer.writeToFile(fileToWrite);

      // Read it back
      pe::UniverseSequence seq = pe::UniverseSequence();
      seq.setDataFileAtomStyle({ atomStyle });
      seq.initializeFromDataSequence({ { fileToWrite } });

      pe::Universe readUniverse = seq.atIndex(0);

      // Verify properties were preserved
      for (size_t i = 1; i <= 5; ++i) {
        pe::Atom prevAtom = universe.getAtom(i);
        pe::Atom readAtom = readUniverse.getAtom(i);

        CHECK(prevAtom == readAtom);
      }

      // Clean up
      std::filesystem::remove(fileToWrite);
    }
  }
}

TEST_CASE("Data-files can be written with velocities",
          "[utils][DataFileWriter][DataFileParser][io]")
{
  std::cout << "Running test \"Data-files can be written with velocities\""
            << std::endl;
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  CHECK(std::filesystem::exists(suspectedPath));

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  std::string largeInputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  universeSeq.initializeFromDataSequence({ { largeInputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  universe.removeAllDihedralAngles();
  universe.removeAllAngles();

  // add velocities
  universe.resampleVelocities(1.0, 3.0);

  // write data file
  pu::DataFileWriter writer = pu::DataFileWriter(universe);
  writer.configIncludeVelocities(true);
  writer.configMoleculeIdxForSwap(false);
  writer.configMoveIntoBox(false);
  writer.configIncludeAngles(false);
  writer.configIncludeDihedralAngles(false);
  writer.configAttemptImageReset(false);
  std::string fileToWrite =
    suspectedPath + "/tmp_data_file_with_velocities.structure.out";

  SECTION("Without re-indexing or moving of images")
  {
    writer.configAtomStyle(pu::AtomStyle::MOLECULAR);
    writer.writeToFile(fileToWrite);

    pe::UniverseSequence seq = pe::UniverseSequence();
    seq.initializeFromDataSequence({ { fileToWrite } });

    pe::Universe readUniverse = seq.atIndex(0);

    CHECK(readUniverse.vertexPropertyExists("vx"));
    CHECK(readUniverse.vertexPropertyExists("vy"));
    CHECK(readUniverse.vertexPropertyExists("vz"));

    CHECK(universe.getBox() == readUniverse.getBox());
    CHECK(universe == readUniverse);
  }

  SECTION("Coordinates into box")
  {
    writer.configAtomStyle(pu::AtomStyle::ANGLE);
    writer.configMoveIntoBox(true);
    writer.writeToFile(fileToWrite);

    pe::UniverseSequence seq = pe::UniverseSequence();
    seq.initializeFromDataSequence({ { fileToWrite } });

    pe::Universe readUniverse = seq.atIndex(0);

    CHECK(readUniverse.vertexPropertyExists("vx"));
    CHECK(readUniverse.vertexPropertyExists("vy"));
    CHECK(readUniverse.vertexPropertyExists("vz"));

    REQUIRE(universe.getBox() == readUniverse.getBox());
    std::vector<pe::Atom> atoms = universe.getAtoms();
    for (size_t i = 0; i < atoms.size(); i++) {
      pe::Atom readAtom = readUniverse.getAtom(atoms[i].getId());
      REQUIRE(atoms[i].getId() == readAtom.getId());
      Eigen::Vector3d prevCoords =
        atoms[i].getUnwrappedCoordinates(universe.getBox());
      Eigen::Vector3d readCoords =
        readAtom.getUnwrappedCoordinates(readUniverse.getBox());

      CHECK(prevCoords.isApprox(readCoords));
    }
  }

  std::filesystem::remove(fileToWrite);
}

TEST_CASE("Atom type FULL can be read/written",
          "[DataFileWriter][DataFileParser][io]")
{
  std::cout << "Running test \"Atom type FULL can be read/written\"";

  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  CHECK(std::filesystem::exists(suspectedPath));

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  std::string largeInputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  universeSeq.initializeFromDataSequence({ { largeInputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  universe.removeAllDihedralAngles();
  universe.removeAllAngles();

  // add velocities
  universe.resampleVelocities(1.0, 3.0);
  // add charge
  universe.setPropertyValue(0, "charge", 1.075);

  // write data file
  pu::DataFileWriter writer = pu::DataFileWriter(universe);
  writer.configAtomStyle(pu::AtomStyle::FULL);
  std::string fileToWrite =
    suspectedPath + "/tmp_data_file_style_full.structure.out";
  writer.writeToFile(fileToWrite);

  // read it again
  pu::DataFileParser parser = pu::DataFileParser();
  parser.read(fileToWrite, pu::AtomStyle::FULL);
  CHECK(parser.getAdditionalAtomData().at("charge").size() ==
        universe.getNrOfAtoms());
  CHECK_THAT(parser.getAdditionalAtomData().at("charge")[0],
             Catch::Matchers::WithinRel(1.075));

  std::filesystem::remove(fileToWrite);
}

TEST_CASE("AveFileReader works", "[AveFileReader][io][utils]")
{
  std::cout << "Running test \"AveFileReader works\"" << std::endl;
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  CHECK(std::filesystem::exists(suspectedPath));
  pu::AveFileReader reader =
    pu::AveFileReader(suspectedPath + "/example_avg_file.out.avg.txt");

  CHECK(reader.getNrOfRows() == 5);
  CHECK(reader.getNrOfColumns() == 3);
  std::vector<std::string> columnNames = reader.getColumnNames();
  CHECK(columnNames.size() == 3);
  CHECK(columnNames[0] == "TimeStep");
  CHECK(columnNames[2] == "else");
  std::vector<std::vector<double>> data = reader.getData();
  CHECK(data[0][0] == 100);
  CHECK(data[2][1] == 6000);
  CHECK(data[0].size() == 5);

  SECTION("Autocorrelation is correct")
  {
    std::vector<size_t> dts = { 1, 2 };
    std::vector<double> results = reader.autocorrelateColumn(2, dts);
    CHECK(results.size() == dts.size());
    // ((5000*6000)+(6000*1000)+(1000*8000)+(8000*9000))/4
    CHECK(results[0] == Catch::Approx(29000000.0));
  }

  SECTION("Autocorrelation works on columns")
  {
    std::vector<double> results =
      reader.autocorrelateColumnDifference(1, 2, { 1, 2 });
    CHECK(results.size() == 2);
  }
}
