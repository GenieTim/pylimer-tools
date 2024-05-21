#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/io/AveFileReader.h"
#include "../../src/pylimer_tools_cpp/io/DataFileParser.h"
#include "../../src/pylimer_tools_cpp/io/DataFileWriter.h"
#include "../../src/pylimer_tools_cpp/io/DumpFileParser.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <filesystem>
#include <fstream>
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
  std::cout << "Running test \"FileParsers can be used\"" << std::endl;
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  CHECK(std::filesystem::exists(suspectedPath));

  SECTION("Reading from dump file works")
  {
    pu::DumpFileParser parser =
      pu::DumpFileParser(suspectedPath + "lammps_dump_small.lammpstrj");
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
      pu::DumpFileParser(suspectedPath + "lammps_dump_small.lammpstrj");
    CHECK(parser.getValuesForAt<double>(0, "BOX BOUNDS", 1).size() == 3);
    CHECK_THROWS(parser.getValuesForAt<double>(0, "NOT EXISTING", 9));

    // test throws
    CHECK_THROWS(pu::DumpFileParser("not-existing-file.out"));

    // test without atoms
    pu::DumpFileParser parser4 = pu::DumpFileParser(
      suspectedPath + "lammps_dump_small_no_atoms.lammpstrj");
    CHECK(parser4.getLength() == 1);
    CHECK(parser.getValuesForAt<double>(0, "BOX BOUNDS", 1).size() == 3);
    CHECK(parser.getValuesForAt<double>(0, "BOX BOUNDS", 1)[0] ==
            4.8545999999999999e+01);
  }

  SECTION("Reading from data files works")
  {
    pu::DataFileParser parser = pu::DataFileParser();
    parser.read(suspectedPath + "lammps_data_file.out");
    CHECK(parser.getNrOfAtoms() == 3000);

    // call copy constructor
    pu::DataFileParser parser2 = pu::DataFileParser();
    parser2 = parser;
    CHECK(parser2.getNrOfAtoms() == 3000);
    // test throws
    CHECK_THROWS(parser.read("not-existing-file.out"));

    // test angles reading
    pu::DataFileParser parser3 = pu::DataFileParser();
    parser3.read(suspectedPath + "lammps_data_file_small_wangles.out");
    CHECK(parser3.getNrOfAngles() == 1);


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
    CHECK_NOTHROW(parser.readNGroups(9, 12));
    CHECK(parser.hasKey("BOX BOUNDS") == true);
    CHECK(parser.hasKey("NO EXISTING") == false);

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
  std::cout << "Running test \"Writers can be used\"" << std::endl;
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  CHECK(std::filesystem::exists(suspectedPath));

  SECTION("Files are read and written")
  {
    // TODO: implement
    pe::UniverseSequence universeSeq = pe::UniverseSequence();
    std::string largeInputFile =
      suspectedPath + "structure/network_100_a_46.structure.out";
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
    std::string fileToWrite = suspectedPath + "tmp_data_file.structure.out";
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

    std::filesystem::remove(fileToWrite);
  }
}

TEST_CASE("AveFileReader works", "[AveFileReader][io][utils]")
{
  std::cout << "Running test \"AveFileReader works\"" << std::endl;
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  CHECK(std::filesystem::exists(suspectedPath));
  pu::AveFileReader reader =
    pu::AveFileReader(suspectedPath + "example_avg_file.out.avg.txt");

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
}
