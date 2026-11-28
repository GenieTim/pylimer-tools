#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/io/AveFileReader.h"
#include "../../src/pylimer_tools_cpp/io/CSVSplitter.h"
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

    // test CRLF line endings and inline comments on group headers
    const std::string crlfDumpPath =
      suspectedPath + "/tmp_crlf_comments_dump.lammpstrj";
    {
      std::ofstream out(crlfDumpPath, std::ios::binary | std::ios::trunc);
      out << "  # leading comment\r\n"
          << "ITEM: TIMESTEP # first step\r\n"
          << "10\r\n"
          << "ITEM: NUMBER OF ATOMS\r\n"
          << "1\r\n"
          << "ITEM: BOX BOUNDS pp pp pp\r\n"
          << "0.0 5.0\r\n"
          << "0.0 6.0\r\n"
          << "0.0 7.0\r\n"
          << "ITEM: ATOMS id type xs ys zs\r\n"
          << "1 2 0.5 0.5 0.5\r\n"
          << "  ITEM: TIMESTEP # second step\r\n"
          << "20\r\n"
          << "ITEM: NUMBER OF ATOMS\r\n"
          << "1\r\n"
          << "ITEM: BOX BOUNDS pp pp pp\r\n"
          << "0.0 5.0\r\n"
          << "0.0 6.0\r\n"
          << "0.0 7.0\r\n"
          << "ITEM: ATOMS id type xs ys zs\r\n"
          << "1 2 0.2 0.5 0.8\r\n";
    }
    {
      pu::DumpFileParser crlfParser(crlfDumpPath);
      REQUIRE(crlfParser.getLength() == 2);
      const std::vector<long int> crlfSteps = crlfParser.readTimeSteps();
      REQUIRE(crlfSteps == std::vector<long int>{ 10, 20 });
      CHECK(crlfParser.getValuesForAt<long int>(1, "TIMESTEP", 0) ==
            std::vector<long int>{ 20 });
      const std::vector<std::vector<pe::Atom>> crlfAtoms =
        crlfParser.readAtoms();
      REQUIRE(crlfAtoms.size() == 2);
      CHECK(crlfAtoms[1][0].getX() == Catch::Approx(1.0));
      CHECK(crlfAtoms[1][0].getY() == Catch::Approx(3.0));
      CHECK(crlfAtoms[1][0].getZ() == Catch::Approx(5.6));
      crlfParser.finish();
    }
    std::filesystem::remove(crlfDumpPath);
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
    std::map<std::string, std::vector<long int>> angles =
      universe.detectAngles();
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
    std::map<std::string, std::vector<long int>> detectedAngles =
      universe.detectDihedralAngles();
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

    std::map<std::string, std::vector<long int>> newBonds =
      readUniverse.getBonds();
    std::map<std::string, std::vector<long int>> oldBonds = universe.getBonds();

    REQUIRE(newBonds["bond_from"].size() == oldBonds["bond_from"].size());
    for (size_t i = 0; i < newBonds["bond_from"].size(); i++) {
      CHECK(newBonds["bond_type"][i] == oldBonds["bond_type"][i]);
      CHECK(std::minmax(newBonds["bond_from"][i], newBonds["bond_to"][i]) ==
            std::minmax(oldBonds["bond_from"][i], oldBonds["bond_to"][i]));
    }

    std::map<std::string, std::vector<long int>> newAngles =
      readUniverse.getAngles();
    std::map<std::string, std::vector<long int>> oldAngles =
      universe.getAngles();

    REQUIRE(newAngles["angle_from"].size() == oldAngles["angle_from"].size());
    for (size_t i = 0; i < newAngles["angle_from"].size(); i++) {
      CHECK(newAngles["angle_type"][i] == oldAngles["angle_type"][i]);
      CHECK(newAngles["angle_from"][i] == oldAngles["angle_from"][i]);
      CHECK(newAngles["angle_via"][i] == oldAngles["angle_via"][i]);
      CHECK(newAngles["angle_to"][i] == oldAngles["angle_to"][i]);
    }

    std::filesystem::remove(fileToWrite);
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
            CHECK_THAT(universe.getAtom(atomId).getExtraData().at(property),
                       Catch::Matchers::WithinAbs(intVal, 1e-9));
          } else {
            // Double properties
            const double doubleVal = realDist(gen);
            universe.setPropertyValue(atomId - 1, property.c_str(), doubleVal);
            CHECK_THAT(universe.getAtom(atomId).getExtraData().at(property),
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

  SECTION("Fresh AveFileReader autocorrelates without prior getData and "
          "handles errors")
  {
    pu::AveFileReader freshReader =
      pu::AveFileReader(suspectedPath + "/example_avg_file.out.avg.txt");
    std::vector<double> r1 = freshReader.autocorrelateColumn(2, { 1 });
    CHECK(r1.size() == 1);

    pu::AveFileReader freshReader2 =
      pu::AveFileReader(suspectedPath + "/example_avg_file.out.avg.txt");
    std::vector<double> r2 =
      freshReader2.autocorrelateColumnDifference(1, 2, { 1 });
    CHECK(r2.size() == 1);

    pu::AveFileReader missingReader("/nonexistent/path/to/avg_file.txt");
    CHECK_THROWS(missingReader.getData());
  }

  SECTION("AveFileReader handles CRLF, blank lines, multi-space/tab headers, "
          "and zero-row files")
  {
    const std::filesystem::path crlfAvePath =
      std::filesystem::temp_directory_path() / "pylimer_test_ave_crlf.txt";
    {
      std::ofstream out(crlfAvePath, std::ios::binary | std::ios::trunc);
      out << "\r\n"
          << "# Time-averaged data\r\n"
          << "#   TimeStep\t  c_temp    v_press  \r\n"
          << "\r\n"
          << "100   1.5   2.5 # inline comment\r\n"
          << "\r\n"
          << "200\t2.0\t3.5\r\n"
          << "   \r\n";
    }
    {
      pu::AveFileReader crlfReader(crlfAvePath.string());
      CHECK(crlfReader.getColumnNames() ==
            std::vector<std::string>{ "TimeStep", "c_temp", "v_press" });
      CHECK(crlfReader.getNrOfColumns() == 3);
      CHECK(crlfReader.getNrOfRows() == 2);
      const std::vector<std::vector<double>>& d = crlfReader.getData();
      REQUIRE(d.size() == 3);
      CHECK(d[0] == std::vector<double>{ 100.0, 200.0 });
      CHECK(d[1] == std::vector<double>{ 1.5, 2.0 });
      CHECK(d[2] == std::vector<double>{ 2.5, 3.5 });
    }
    std::filesystem::remove(crlfAvePath);

    const std::filesystem::path zeroRowAvePath =
      std::filesystem::temp_directory_path() /
      "pylimer_test_ave_empty_rows.txt";
    {
      std::ofstream out(zeroRowAvePath, std::ios::binary | std::ios::trunc);
      out << "# Header\r\n# Step Val\r\n\r\n";
    }
    {
      pu::AveFileReader zeroReader(zeroRowAvePath.string());
      CHECK(zeroReader.getNrOfRows() == 0);
      CHECK(zeroReader.getData()[0].empty());
    }
    std::filesystem::remove(zeroRowAvePath);
  }
}

TEST_CASE("CSVSplitter splits without leaking temp files",
          "[CSVSplitter][io][utils]")
{
  const std::filesystem::path tempDir = std::filesystem::temp_directory_path();
  auto countTempCsvFiles = [&tempDir]() -> size_t {
    size_t count = 0;
    for (const std::filesystem::directory_entry& entry :
         std::filesystem::directory_iterator(tempDir)) {
      const std::string name = entry.path().filename().string();
      if (name.starts_with("pylimer_csv_") && name.ends_with(".tmp")) {
        ++count;
      }
    }
    return count;
  };

  const size_t initialTempCount = countTempCsvFiles();

  // Empty file should not create any temporary files
  const std::filesystem::path emptyInput =
    tempDir / "pylimer_test_csv_empty.txt";
  {
    std::ofstream out(emptyInput);
  }
  const std::vector<std::string> emptyResults =
    pu::splitCSV(emptyInput.string(), " ");
  CHECK(emptyResults.empty());
  CHECK(countTempCsvFiles() == initialTempCount);
  std::filesystem::remove(emptyInput);

  // Multi-section file should create only the returned temporary files
  const std::filesystem::path multiInput =
    tempDir / "pylimer_test_csv_multi.txt";
  {
    std::ofstream out(multiInput);
    out << "Step a b\n";
    out << "1 2 3\n";
    out << "Step a b c\n";
    out << "1 2 3 4\n";
  }
  const std::vector<std::string> results =
    pu::splitCSV(multiInput.string(), " ");
  CHECK(results.size() == 2);
  CHECK(countTempCsvFiles() == initialTempCount + results.size());

  for (const std::string& tmpPath : results) {
    CHECK(std::filesystem::exists(tmpPath));
    std::filesystem::remove(tmpPath);
  }
  std::filesystem::remove(multiInput);
}

namespace {
struct TempFileGuard
{
  std::filesystem::path path;
  explicit TempFileGuard(const std::string& filename)
    : path(std::filesystem::temp_directory_path() / filename)
  {
  }
  ~TempFileGuard()
  {
    std::error_code ec;
    std::filesystem::remove(path, ec);
  }
  void writeBinary(const std::string& content) const
  {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    REQUIRE(out.is_open());
    out.write(content.data(), static_cast<std::streamsize>(content.size()));
    REQUIRE(out.good());
  }
};
} // namespace

TEST_CASE("DumpFileParser handles boundary cases",
          "[utils][DumpFileParser][io]")
{
  SECTION("No trailing newline at EOF and zero-atom timestep")
  {
    TempFileGuard tmp("pylimer_test_dump_no_eof_nl.lammpstrj");
    tmp.writeBinary("ITEM: TIMESTEP\r\n"
                    "0\r\n"
                    "ITEM: NUMBER OF ATOMS\r\n"
                    "0\r\n"
                    "ITEM: BOX BOUNDS pp pp pp\r\n"
                    "-1.5 2.5\r\n"
                    "0.0 4.0\r\n"
                    "1.0e-1 3.1e+0\r\n"
                    "ITEM: ATOMS id type x y z\r\n"
                    "ITEM: TIMESTEP\r\n"
                    "5\r\n"
                    "ITEM: NUMBER OF ATOMS\r\n"
                    "1\r\n"
                    "ITEM: BOX BOUNDS pp pp pp\r\n"
                    "0.0 10.0\r\n"
                    "0.0 10.0\r\n"
                    "0.0 10.0\r\n"
                    "ITEM: ATOMS id type x y z\r\n"
                    "7 3 1.25 2.5 3.75");

    pu::DumpFileParser parser(tmp.path.string());
    REQUIRE(parser.getLength() == 2);
    CHECK(parser.readTimeSteps() == std::vector<long int>{ 0, 5 });

    const std::vector<pe::Box> boxes = parser.readBoxes();
    REQUIRE(boxes.size() == 2);
    CHECK(boxes[0].getLx() == Catch::Approx(4.0));
    CHECK(boxes[0].getLy() == Catch::Approx(4.0));
    CHECK(boxes[0].getLz() == Catch::Approx(3.0));

    const std::vector<std::vector<pe::Atom>> atoms = parser.readAtoms();
    REQUIRE(atoms.size() == 2);
    CHECK(atoms[0].empty());
    REQUIRE(atoms[1].size() == 1);
    CHECK(atoms[1][0].getId() == 7);
    CHECK(atoms[1][0].getType() == 3);
    CHECK(atoms[1][0].getX() == Catch::Approx(1.25));
    CHECK(atoms[1][0].getY() == Catch::Approx(2.5));
    CHECK(atoms[1][0].getZ() == Catch::Approx(3.75));
  }

  SECTION("CRLF and LF produce identical results across coordinate formats")
  {
    const std::string lfContent = "ITEM: TIMESTEP\n"
                                  "100\n"
                                  "ITEM: NUMBER OF ATOMS\n"
                                  "2\n"
                                  "ITEM: BOX BOUNDS pp pp pp\n"
                                  "0.0 10.0\n"
                                  "0.0 20.0\n"
                                  "0.0 30.0\n"
                                  "ITEM: ATOMS type z y x id\n"
                                  "2 3.0 2.0 1.0 10\n"
                                  "1 6.0 4.0 2.0 20\n"
                                  "ITEM: TIMESTEP\n"
                                  "200\n"
                                  "ITEM: NUMBER OF ATOMS\n"
                                  "1\n"
                                  "ITEM: BOX BOUNDS pp pp pp\n"
                                  "0.0 10.0\n"
                                  "0.0 20.0\n"
                                  "0.0 30.0\n"
                                  "ITEM: ATOMS id type xsu ysu zsu vx\n"
                                  "10 2 0.5 0.25 0.1 1.5\n";

    std::string crlfContent;
    crlfContent.reserve(lfContent.size() + 32);
    for (const char c : lfContent) {
      if (c == '\n') {
        crlfContent.push_back('\r');
      }
      crlfContent.push_back(c);
    }

    TempFileGuard lfTmp("pylimer_test_dump_lf.lammpstrj");
    TempFileGuard crlfTmp("pylimer_test_dump_crlf.lammpstrj");
    lfTmp.writeBinary(lfContent);
    crlfTmp.writeBinary(crlfContent);

    pu::DumpFileParser lfParser(lfTmp.path.string());
    pu::DumpFileParser crlfParser(crlfTmp.path.string());

    REQUIRE(lfParser.getLength() == 2);
    REQUIRE(crlfParser.getLength() == 2);

    const pu::ReadDumpFileSections lfSections = lfParser.readDumpFileSections(
      pu::ReadableDumpFileSections::TIMESTEP |
      pu::ReadableDumpFileSections::BOX | pu::ReadableDumpFileSections::ATOM |
      pu::ReadableDumpFileSections::EXTRA_ATOM);
    const pu::ReadDumpFileSections crlfSections =
      crlfParser.readDumpFileSections(pu::ReadableDumpFileSections::TIMESTEP |
                                      pu::ReadableDumpFileSections::BOX |
                                      pu::ReadableDumpFileSections::ATOM |
                                      pu::ReadableDumpFileSections::EXTRA_ATOM);

    CHECK(lfSections.timesteps == crlfSections.timesteps);
    REQUIRE(lfSections.atoms.size() == 2);
    REQUIRE(crlfSections.atoms.size() == 2);
    CHECK(lfSections.atoms[0][0] == crlfSections.atoms[0][0]);
    CHECK(lfSections.atoms[0][1] == crlfSections.atoms[0][1]);
    CHECK(lfSections.atoms[1][0] == crlfSections.atoms[1][0]);
    CHECK(crlfSections.atoms[1][0].getX() == Catch::Approx(5.0));
    CHECK(crlfSections.atoms[1][0].getY() == Catch::Approx(5.0));
    CHECK(crlfSections.atoms[1][0].getZ() == Catch::Approx(3.0));
    CHECK(lfSections.extraAtomsData[1].at("vx") ==
          crlfSections.extraAtomsData[1].at("vx"));
    CHECK(crlfSections.extraAtomsData[1].at("vx") ==
          std::vector<double>{ 1.5 });
  }

  SECTION("Chunk boundary (>256 KiB) with CRLF and random group access")
  {
    TempFileGuard tmp("pylimer_test_dump_large_crlf.lammpstrj");
    std::string content;
    content.reserve(300 * 1024);
    // Step 0: enough atoms to exceed 256 KiB scan buffer
    const int nAtoms = 12000;
    content += "ITEM: TIMESTEP\r\n0\r\n";
    content += "ITEM: NUMBER OF ATOMS\r\n" + std::to_string(nAtoms) + "\r\n";
    content +=
      "ITEM: BOX BOUNDS pp pp pp\r\n0.0 10.0\r\n0.0 10.0\r\n0.0 10.0\r\n";
    content += "ITEM: ATOMS id type x y z ix iy iz\r\n";
    for (int i = 1; i <= nAtoms; ++i) {
      content += std::to_string(i) + " 1 1.0 2.0 3.0 0 0 0\r\n";
    }
    // Step 1 after 256 KiB boundary
    content += "ITEM: TIMESTEP\r\n42\r\n";
    content += "ITEM: NUMBER OF ATOMS\r\n1\r\n";
    content +=
      "ITEM: BOX BOUNDS pp pp pp\r\n0.0 10.0\r\n0.0 10.0\r\n0.0 10.0\r\n";
    content += "ITEM: ATOMS id type x y z ix iy iz\r\n";
    content += "9999 2 4.0 5.0 6.0 1 -1 2\r\n";
    REQUIRE(content.size() > 260 * 1024);
    tmp.writeBinary(content);

    pu::DumpFileParser parser(tmp.path.string());
    REQUIRE(parser.getLength() == 2);
    CHECK(parser.readTimeSteps() == std::vector<long int>{ 0, 42 });

    // Random access in reverse order + forgetAt
    parser.readGroupByIdx(1);
    CHECK(parser.getValuesForAt<long int>(1, "ATOMS", "id") ==
          std::vector<long int>{ 9999 });
    CHECK(parser.getValuesForAt<int>(1, "ATOMS", "ix") ==
          std::vector<int>{ 1 });
    parser.forgetAt(1);
    CHECK(parser.getValuesForAt<long int>(0, "TIMESTEP", 0) ==
          std::vector<long int>{ 0 });
    CHECK(parser.getValuesForAt<long int>(1, "TIMESTEP", 0) ==
          std::vector<long int>{ 42 });
  }

  SECTION("Truncated and closed-file error handling")
  {
    TempFileGuard truncatedTimestep("pylimer_test_dump_trunc_ts.lammpstrj");
    truncatedTimestep.writeBinary("ITEM: TIMESTEP\r\n");
    pu::DumpFileParser p1(truncatedTimestep.path.string());
    CHECK_THROWS_AS(p1.readTimeSteps(), std::runtime_error);

    TempFileGuard truncatedBox("pylimer_test_dump_trunc_box.lammpstrj");
    truncatedBox.writeBinary(
      "ITEM: TIMESTEP\n0\nITEM: BOX BOUNDS pp pp pp\n0 1\n");
    pu::DumpFileParser p2(truncatedBox.path.string());
    CHECK_THROWS_AS(p2.readBoxes(), std::runtime_error);

    TempFileGuard truncatedAtoms("pylimer_test_dump_trunc_atoms.lammpstrj");
    truncatedAtoms.writeBinary("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
                               "ITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\n"
                               "ITEM: ATOMS id type x y z\n1 1 0.1 0.2 0.3\n");
    pu::DumpFileParser p3(truncatedAtoms.path.string());
    CHECK_THROWS_AS(p3.readAtoms(), std::runtime_error);

    p3.finish();
    CHECK_THROWS_AS(p3.readGroupByIdx(0), std::runtime_error);

    TempFileGuard emptyDump("pylimer_test_dump_empty.lammpstrj");
    emptyDump.writeBinary("   # comment only\r\n\r\n");
    pu::DumpFileParser pEmpty(emptyDump.path.string());
    CHECK(pEmpty.getLength() == 0);
    CHECK(pEmpty.readTimeSteps().empty());
    CHECK_NOTHROW(pEmpty.read());

    TempFileGuard malformedAtom("pylimer_test_dump_bad_atom.lammpstrj");
    malformedAtom.writeBinary("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
                              "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
                              "ITEM: ATOMS id type x y z ix iy iz\n"
                              "1 1 0.5 0.5\n");
    pu::DumpFileParser pBadAtom(malformedAtom.path.string());
    CHECK_THROWS_AS(pBadAtom.readAtoms(), std::runtime_error);
  }

  SECTION("Scaled coordinates respect non-zero lower box bounds and 64-bit IDs")
  {
    TempFileGuard tmp("pylimer_test_dump_shifted_box.lammpstrj");
    tmp.writeBinary("ITEM: TIMESTEP\n10\n"
                    "ITEM: NUMBER OF ATOMS\n1\n"
                    "ITEM: BOX BOUNDS pp pp pp\n"
                    "-5.0 5.0\n"
                    "2.0 12.0\n"
                    "-10.0 0.0\n"
                    "ITEM: ATOMS id type xs ys zs\n"
                    "5000000000 2 0.5 0.25 0.8\n");

    pu::DumpFileParser parser(tmp.path.string());
    const std::vector<std::vector<pe::Atom>> atoms = parser.readAtoms();
    REQUIRE(atoms.size() == 1);
    REQUIRE(atoms[0].size() == 1);
    CHECK(atoms[0][0].getId() == 5000000000L);
    CHECK(atoms[0][0].getType() == 2);
    CHECK(atoms[0][0].getX() == Catch::Approx(0.0));
    CHECK(atoms[0][0].getY() == Catch::Approx(4.5));
    CHECK(atoms[0][0].getZ() == Catch::Approx(-2.0));
  }
}

TEST_CASE("DataFileParser handles boundary cases",
          "[utils][DataFileParser][io]")
{
  SECTION("CRLF vs LF equivalence with comments, tabs, and all sections")
  {
    const std::string lfContent = "LAMMPS data file comment\n"
                                  "\n"
                                  "3 atoms # total atoms\n"
                                  "2 bonds\n"
                                  "1 angles\n"
                                  "0 dihedrals\n"
                                  "2 atom types\n"
                                  "1 bond types\n"
                                  "1 angle types\n"
                                  "\n"
                                  "-5.0 5.0 xlo xhi # x bounds\n"
                                  "0.0 1.0e1 ylo yhi\n"
                                  "-2.5 7.5 zlo zhi\n"
                                  "\n"
                                  "Masses # bead masses\n"
                                  "\n"
                                  "1 1.0 # type 1\n"
                                  "2\t2.5 # type 2\n"
                                  "\n"
                                  "Atoms # bond\n"
                                  "\n"
                                  "1 1 1 0.5 1.5 2.5 0 1 -1\n"
                                  "2\t1\t2\t1.5\t2.5\t3.5\t1\t0\t0\n"
                                  "3 2 1 -1.0 0.0 4.0 0 0 1\n"
                                  "\n"
                                  "Velocities # out of atom-id order\n"
                                  "\n"
                                  "3 0.3 0.6 0.9\n"
                                  "1 0.1 0.2 0.3\n"
                                  "2 0.2 0.4 0.6\n"
                                  "\n"
                                  "Bonds\n"
                                  "\n"
                                  "1 1 1 2\n"
                                  "2 1 2 3\n"
                                  "\n"
                                  "Angles\n"
                                  "\n"
                                  "1 1 1 2 3\n";

    std::string crlfContent;
    crlfContent.reserve(lfContent.size() + 64);
    for (const char c : lfContent) {
      if (c == '\n') {
        crlfContent.push_back('\r');
      }
      crlfContent.push_back(c);
    }

    TempFileGuard lfTmp("pylimer_test_data_lf.out");
    TempFileGuard crlfTmp("pylimer_test_data_crlf.out");
    lfTmp.writeBinary(lfContent);
    crlfTmp.writeBinary(crlfContent);

    pu::DataFileParser lfParser;
    pu::DataFileParser crlfParser;
    lfParser.read(lfTmp.path.string());
    crlfParser.read(crlfTmp.path.string());

    CHECK(crlfParser.getNrOfAtoms() == 3);
    CHECK(crlfParser.getNrOfBonds() == 2);
    CHECK(crlfParser.getNrOfAngles() == 1);
    CHECK(crlfParser.getNrOfDihedralAngles() == 0);
    CHECK(crlfParser.getLowX() == Catch::Approx(-5.0));
    CHECK(crlfParser.getHighX() == Catch::Approx(5.0));
    CHECK(crlfParser.getLx() == Catch::Approx(10.0));
    CHECK(crlfParser.getLy() == Catch::Approx(10.0));
    CHECK(crlfParser.getLz() == Catch::Approx(10.0));

    CHECK(crlfParser.getAtomIds() == lfParser.getAtomIds());
    CHECK(crlfParser.getMoleculeIds() == lfParser.getMoleculeIds());
    CHECK(crlfParser.getAtomTypes() == lfParser.getAtomTypes());
    CHECK(crlfParser.getAtomX() == lfParser.getAtomX());
    CHECK(crlfParser.getAtomY() == lfParser.getAtomY());
    CHECK(crlfParser.getAtomZ() == lfParser.getAtomZ());
    CHECK(crlfParser.getAtomNx() == lfParser.getAtomNx());
    CHECK(crlfParser.getAtomNy() == lfParser.getAtomNy());
    CHECK(crlfParser.getAtomNz() == lfParser.getAtomNz());
    CHECK(crlfParser.getMasses() == lfParser.getMasses());
    CHECK(crlfParser.getBondFrom() == lfParser.getBondFrom());
    CHECK(crlfParser.getBondTo() == lfParser.getBondTo());
    CHECK(crlfParser.getAngleFrom() == lfParser.getAngleFrom());
    CHECK(crlfParser.getAngleVia() == lfParser.getAngleVia());
    CHECK(crlfParser.getAngleTo() == lfParser.getAngleTo());
    // Velocities must be reordered to match atomIds [1, 2, 3]
    CHECK(crlfParser.getAdditionalAtomData().at("vx") ==
          std::vector<double>{ 0.1, 0.2, 0.3 });
    CHECK(crlfParser.getAdditionalAtomData().at("vy") ==
          std::vector<double>{ 0.2, 0.4, 0.6 });
    CHECK(crlfParser.getAdditionalAtomData().at("vz") ==
          std::vector<double>{ 0.3, 0.6, 0.9 });
  }

  SECTION("No trailing newline at EOF, zero bonds, and explicit atom style")
  {
    TempFileGuard tmp("pylimer_test_data_no_eof_nl.out");
    tmp.writeBinary("Header\r\n"
                    "2 atoms\r\n"
                    "1 atom types\r\n"
                    "0.0 5.0 xlo xhi\r\n"
                    "0.0 5.0 ylo yhi\r\n"
                    "0.0 5.0 zlo zhi\r\n"
                    "\r\n"
                    "Masses\r\n"
                    "\r\n"
                    "1 1.0\r\n"
                    "\r\n"
                    "Atoms\r\n"
                    "\r\n"
                    "1 1 0.1 0.2 0.3 0 0 0\r\n"
                    "2 1 1.1 1.2 1.3 0 0 0");

    pu::DataFileParser parser;
    parser.read(tmp.path.string(), pu::AtomStyle::ATOMIC);
    CHECK(parser.getNrOfAtoms() == 2);
    CHECK(parser.getNrOfBonds() == 0);
    CHECK(parser.getAtomIds() == std::vector<long int>{ 1, 2 });
    CHECK(parser.getAtomX() == std::vector<double>{ 0.1, 1.1 });
  }

  SECTION("Malformed and truncated data files throw descriptive errors")
  {
    TempFileGuard badMass("pylimer_test_data_bad_mass.out");
    badMass.writeBinary("Header\n1 atoms\n1 atom types\n"
                        "0 1 xlo xhi\n0 1 ylo yhi\n0 1 zlo zhi\n\n"
                        "Masses\n\n1 1.0 2.0\n\nAtoms\n\n1 1 1 0 0 0 0 0 0\n");
    pu::DataFileParser p1;
    CHECK_THROWS_AS(p1.read(badMass.path.string()), std::runtime_error);

    TempFileGuard truncAtoms("pylimer_test_data_trunc_atoms.out");
    truncAtoms.writeBinary(
      "Header\n2 atoms\n1 atom types\n"
      "0 1 xlo xhi\n0 1 ylo yhi\n0 1 zlo zhi\n\n"
      "Masses\n\n1 1.0\n\nAtoms # atomic\n\n1 1 0 0 0 0 0 0\n");
    pu::DataFileParser p2;
    CHECK_THROWS_AS(p2.read(truncAtoms.path.string()), std::runtime_error);

    TempFileGuard truncBonds("pylimer_test_data_trunc_bonds.out");
    truncBonds.writeBinary(
      "Header\n2 atoms\n2 bonds\n1 atom types\n1 bond types\n"
      "0 1 xlo xhi\n0 1 ylo yhi\n0 1 zlo zhi\n\n"
      "Masses\n\n1 1.0\n\nAtoms # atomic\n\n"
      "1 1 0 0 0 0 0 0\n2 1 0 0 0 0 0 0\n\nBonds\n\n1 1 1 2\n");
    pu::DataFileParser p3;
    CHECK_THROWS_AS(p3.read(truncBonds.path.string()), std::runtime_error);
  }

  SECTION("Auto-detected Atoms # style, optional Masses, section order, and "
          "parser reuse")
  {
    TempFileGuard autoFullTmp("pylimer_test_data_autofull.out");
    autoFullTmp.writeBinary(
      "Title mentioning 999 atoms and Masses keyword\n"
      "2 atoms\n"
      "0 bonds\n"
      "1 angles\n"
      "1 atom types\n"
      "1 angle types\n"
      "-5.0 5.0 xlo xhi\n"
      "-5.0 5.0 ylo yhi\n"
      "-5.0 5.0 zlo zhi\n\n"
      "Pair Coeffs # lj/cut\n\n"
      "1 1.0 1.0\n\n"
      "Velocities\n\n"
      "2 0.4 0.5 0.6\n"
      "1 0.1 0.2 0.3\n\n"
      "# Bonds in a comment must not trigger bond parsing\n"
      "Bonds\n\n"
      "Atoms # full\n\n"
      "1 10 1 -0.5 1.0 2.0 3.0\n"
      "2 10 1 +0.5 -1.0 -2.0 -3.0\n\n"
      "Angles\n\n"
      "7 1 1 2 1\n");

    pu::DataFileParser parser;
    // Default atomStyle is ANGLE, so # full auto-detection must kick in
    parser.read(autoFullTmp.path.string());
    REQUIRE(parser.getNrOfAtoms() == 2);
    CHECK(parser.getNrOfBonds() == 0);
    CHECK(parser.getBondIds().empty());
    REQUIRE(parser.getNrOfAngles() == 1);
    CHECK(parser.getAngleIds() == std::vector<long int>{ 7 });
    CHECK(parser.getMoleculeIds() == std::vector<int>{ 10, 10 });
    CHECK(parser.getAdditionalAtomData().at("charge") ==
          std::vector<double>{ -0.5, 0.5 });
    CHECK(parser.getAtomNx() == std::vector<int>{ 0, 0 });
    CHECK(parser.getAdditionalAtomData().at("vx") ==
          std::vector<double>{ 0.1, 0.4 });

    // Reusing the same parser instance must reset all previous state
    TempFileGuard hybridTmp("pylimer_test_data_hybrid.out");
    hybridTmp.writeBinary("Hybrid data file\n"
                          "1 atoms\n"
                          "1 atom types\n"
                          "0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\n"
                          "Masses\n\n1 1.0\n\n"
                          "Atoms # hybrid bond edpd\n\n"
                          "5 1 1.0 2.0 3.0 7 300.0 1.25\n");

    parser.read(hybridTmp.path.string());
    REQUIRE(parser.getNrOfAtoms() == 1);
    CHECK(parser.getNrOfAngles() == 0);
    CHECK(parser.getAngleIds().empty());
    CHECK(parser.getAtomIds() == std::vector<long int>{ 5 });
    CHECK(parser.getMoleculeIds() == std::vector<int>{ 7 });
    CHECK(parser.getAtomNx() == std::vector<int>{ 0 });
    CHECK(parser.getAdditionalAtomData().find("charge") ==
          parser.getAdditionalAtomData().end());
    CHECK(parser.getAdditionalAtomData().at("edpd_temp") ==
          std::vector<double>{ 300.0 });
    CHECK(parser.getAdditionalAtomData().at("edpd") ==
          std::vector<double>{ 1.25 });
  }
}

TEST_CASE("IO readers benchmark", "[benchmark][io_benchmark][.benchmark]")
{
  const size_t nDataAtoms = 50000;
  TempFileGuard dataTmp("pylimer_bench_data.out");
  {
    std::string buf;
    buf.reserve(nDataAtoms * 120);
    buf += "LAMMPS data file\n\n";
    buf += std::to_string(nDataAtoms) + " atoms\n";
    buf += std::to_string(nDataAtoms - 1) + " bonds\n";
    buf += "1 atom types\n1 bond types\n";
    buf += "0.0 100.0 xlo xhi\n0.0 100.0 ylo yhi\n0.0 100.0 zlo zhi\n\n";
    buf += "Masses\n\n1 1.0\n\nAtoms # full\n\n";
    for (size_t i = 1; i <= nDataAtoms; ++i) {
      buf += std::to_string(i) + " 1 1 0.1 1.25 2.5 3.75 0 0 0\n";
    }
    buf += "\nBonds\n\n";
    for (size_t i = 1; i < nDataAtoms; ++i) {
      buf += std::to_string(i) + " 1 " + std::to_string(i) + " " +
             std::to_string(i + 1) + "\n";
    }
    buf += "\nVelocities\n\n";
    for (size_t i = 1; i <= nDataAtoms; ++i) {
      buf += std::to_string(i) + " 0.01 -0.02 0.03\n";
    }
    dataTmp.writeBinary(buf);
  }

  const size_t nFrames = 5;
  const size_t nDumpAtoms = 20000;
  TempFileGuard dumpStdTmp("pylimer_bench_dump_std.lammpstrj");
  TempFileGuard dumpExtraTmp("pylimer_bench_dump_extra.lammpstrj");
  {
    std::string bufStd;
    std::string bufExtra;
    bufStd.reserve(nFrames * nDumpAtoms * 45);
    bufExtra.reserve(nFrames * nDumpAtoms * 50);
    for (size_t f = 0; f < nFrames; ++f) {
      const std::string hdr = "ITEM: TIMESTEP\n" + std::to_string(f * 1000) +
                              "\nITEM: NUMBER OF ATOMS\n" +
                              std::to_string(nDumpAtoms) +
                              "\nITEM: BOX BOUNDS pp pp pp\n"
                              "0.0 100.0\n0.0 100.0\n0.0 100.0\n";
      bufStd += hdr + "ITEM: ATOMS id type x y z ix iy iz\n";
      bufExtra += hdr + "ITEM: ATOMS id type xu yu zu vx vy vz\n";
      for (size_t i = 1; i <= nDumpAtoms; ++i) {
        bufStd += std::to_string(i) + " 1 1.25 2.5 3.75 0 1 -1\n";
        bufExtra += std::to_string(i) + " 1 1.25 2.5 3.75 0.1 -0.2 0.3\n";
      }
    }
    dumpStdTmp.writeBinary(bufStd);
    dumpExtraTmp.writeBinary(bufExtra);
  }

  TempFileGuard dumpSingleTmp("pylimer_bench_dump_single.lammpstrj");
  {
    std::string bufSingle;
    bufSingle.reserve(nDumpAtoms * 45);
    bufSingle += "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n" +
                 std::to_string(nDumpAtoms) +
                 "\nITEM: BOX BOUNDS pp pp pp\n"
                 "0.0 100.0\n0.0 100.0\n0.0 100.0\n"
                 "ITEM: ATOMS id type x y z ix iy iz\n";
    for (size_t i = 1; i <= nDumpAtoms; ++i) {
      bufSingle += std::to_string(i) + " 1 1.25 2.5 3.75 0 1 -1\n";
    }
    dumpSingleTmp.writeBinary(bufSingle);
  }

  const size_t nAveRows = 20000;
  TempFileGuard aveTmp("pylimer_bench_ave.txt");
  {
    std::string bufAve;
    bufAve.reserve(nAveRows * 45);
    bufAve += "# Time-averaged data\n# TimeStep c1 c2 c3 c4 c5\n";
    for (size_t i = 0; i < nAveRows; ++i) {
      bufAve += std::to_string(i * 10) + " 1.1 2.2 3.3 4.4 5.5\n";
    }
    aveTmp.writeBinary(bufAve);
  }

  BENCHMARK("DataFileParser::read (50k atoms full + bonds + velocities)")
  {
    pu::DataFileParser p;
    p.read(dataTmp.path.string(), pu::AtomStyle::FULL);
    return p.getNrOfAtoms();
  };

  BENCHMARK("DumpFileParser::readAtoms (5x20k standard format)")
  {
    pu::DumpFileParser p(dumpStdTmp.path.string());
    return p.readAtoms().size();
  };

  BENCHMARK("DumpFileParser::readDumpFileSections (5x20k extra cols)")
  {
    pu::DumpFileParser p(dumpExtraTmp.path.string());
    return p
      .readDumpFileSections(pu::ReadableDumpFileSections::TIMESTEP |
                            pu::ReadableDumpFileSections::BOX |
                            pu::ReadableDumpFileSections::ATOM |
                            pu::ReadableDumpFileSections::EXTRA_ATOM)
      .atoms.size();
  };

  BENCHMARK("DumpFileParser::read + getValuesForAt (20k atoms)")
  {
    pu::DumpFileParser p(dumpSingleTmp.path.string());
    p.read();
    return p.getValuesForAt<double>(0, "ATOMS", "x").size();
  };

  BENCHMARK("AveFileReader::getData (20k rows x 6 cols)")
  {
    pu::AveFileReader r(aveTmp.path.string());
    return r.getData()[0].size();
  };
}
