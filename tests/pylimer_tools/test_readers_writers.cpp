#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/io/DataFileParser.h"
#include "../../src/pylimer_tools_cpp/io/DataFileWriter.h"
#include "../../src/pylimer_tools_cpp/io/DumpFileParser.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
extern "C"
{
#include <igraph/igraph.h>
}

namespace pu = pylimer_tools::utils;
namespace pe = pylimer_tools::entities;

TEST_CASE("FileParsers can be used", "[utils][DumpFileParser][DataFileParser]")
{
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));

  SECTION("Reading from dump file works")
  {
    pu::DumpFileParser parser =
      pu::DumpFileParser(suspectedPath + "lammps_dump_small.lammpstrj");
    REQUIRE_THROWS(parser.hasKey("BOX BOUNDS"));
    REQUIRE(parser.getLength() == 1);
    REQUIRE_NOTHROW(parser.read());
    REQUIRE(parser.getLength() == 1);
    REQUIRE_THROWS(parser.readNGroups(1, 10));
    REQUIRE_THROWS(parser.readNGroups(10, 1));
    REQUIRE_THROWS(parser.getValuesForAt<int>(0, "BOX BOUNDS", "notexisting"));
    // call copy constructor
    pu::DumpFileParser parser2 = parser;
    REQUIRE(parser2.getLength() == 1);

    // test other reading capabilities
    pu::DumpFileParser parser3 =
      pu::DumpFileParser(suspectedPath + "lammps_dump_small.lammpstrj");
    REQUIRE(parser.getValuesForAt<double>(0, "BOX BOUNDS", 1).size() == 3);
    REQUIRE_THROWS(parser.getValuesForAt<double>(0, "NOT EXISTING", 9));

    // test throws
    REQUIRE_THROWS(pu::DumpFileParser("not-existing-file.out"));

    // test without atoms
    pu::DumpFileParser parser4 = pu::DumpFileParser(
      suspectedPath + "lammps_dump_small_no_atoms.lammpstrj");
    REQUIRE(parser4.getLength() == 1);
    REQUIRE(parser.getValuesForAt<double>(0, "BOX BOUNDS", 1).size() == 3);
    REQUIRE(parser.getValuesForAt<double>(0, "BOX BOUNDS", 1)[0] ==
            4.8545999999999999e+01);
  }

  SECTION("Reading from data files works")
  {
    pu::DataFileParser parser = pu::DataFileParser();
    parser.read(suspectedPath + "lammps_data_file.out");
    REQUIRE(parser.getNrOfAtoms() == 3000);
    // call copy constructor
    pu::DataFileParser parser2 = pu::DataFileParser();
    parser2 = parser;
    REQUIRE(parser2.getNrOfAtoms() == 3000);
    // test throws
    REQUIRE_THROWS(parser.read("not-existing-file.out"));

    // test angles reading
    pu::DataFileParser parser3 = pu::DataFileParser();
    parser3.read(suspectedPath + "lammps_data_file_small_wangles.out");
    REQUIRE(parser3.getNrOfAngles() == 1);

    // BENCHMARK("DataFileParserOld")
    // {
    //   pu::DataFileParser parser4 = pu::DataFileParser();
    //   parser4.read(suspectedPath + "big_dump_file_data.out");
    //   return parser4.getNrOfAngles();
    // };

    // BENCHMARK("DataFileParserNew")
    // {
    //   pu::DataFileParser2 parser5 = pu::DataFileParser2();
    //   parser5.read(suspectedPath + "big_dump_file_data.out");
    //   return parser5.getNrOfAngles();
    // };
  }

  SECTION("Reading large files is sensibly fast")
  {
    pu::DumpFileParser parser =
      pu::DumpFileParser(suspectedPath + "big_dump_file.lammpstrj");
    // pre-read multiple
    REQUIRE_NOTHROW(parser.readNGroups(9, 12));
    REQUIRE(parser.hasKey("BOX BOUNDS") == true);
    REQUIRE(parser.hasKey("NO EXISTING") == false);

    // pu::DumpFileParser2 parser5 =
    //   pu::DumpFileParser2(suspectedPath + "big_dump_file.lammpstrj");
    // CHECK(parser5.getLength() == parser.getLength());

    // BENCHMARK("DumpFileParserOld")
    // {
    //   pu::DumpFileParser parser4 =
    //     pu::DumpFileParser(suspectedPath + "big_dump_file.lammpstrj");
    //   return parser4.getLength();
    // };

    // BENCHMARK("DumpFileParserNew")
    // {
    //   pu::DumpFileParser2 parser5 =
    //     pu::DumpFileParser2(suspectedPath + "big_dump_file.lammpstrj");
    //   return parser5.getLength();
    // };
  }
};

TEST_CASE("Writers can be used", "[utils][DataFileWriter][DataFileParser]")
{
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));

  SECTION("Files are read and written")
  {
    // TODO: implement
    pe::UniverseSequence universeSeq = pe::UniverseSequence();
    std::string largeInputFile =
      suspectedPath + "network_83_a_100.structure.out";
    universeSeq.initializeFromDataSequence({ { largeInputFile } });
    pe::Universe universe = universeSeq.atIndex(0);

    // add angles
    auto angles = universe.detectAngles();
    std::vector<int> angleTypes;
    angleTypes.reserve(angles["angle_from"].size());
    for (size_t i = 0; i < angles["angle_from"].size(); i++) {
      angleTypes.push_back(1);
    }
    universe.addAngles(angles["angle_from"],
                       angles["angle_via"],
                       angles["angle_to"],
                       angleTypes);
    REQUIRE(universe.getNrOfAngles() > 0);

    // add dihedrals
    auto detectedAngles = universe.detectDihedralAngles();
    REQUIRE(detectedAngles["dihedral_angle_from"].size() == 8);
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

    // write data file
    pu::DataFileWriter writer = pu::DataFileWriter(universe);
    writer.configIncludeAngles(true);
    writer.configIncludeDihedralAngles(true);
    std::string fileToWrite = suspectedPath + "tmp_data_file.structure.out";
    writer.writeToFile(fileToWrite);

    pe::UniverseSequence seq = pe::UniverseSequence();
    seq.initializeFromDataSequence({ { fileToWrite } });
    pe::Universe readUniverse = seq.atIndex(0);

    REQUIRE(universe.getNrOfAtoms() == readUniverse.getNrOfAtoms());
    REQUIRE(universe.getNrOfBonds() == readUniverse.getNrOfBonds());
    REQUIRE(universe.getNrOfAngles() == readUniverse.getNrOfAngles());
    REQUIRE(universe.getNrOfAngles() == readUniverse.getNrOfDihedralAngles());

    std::remove(fileToWrite);
  }
}