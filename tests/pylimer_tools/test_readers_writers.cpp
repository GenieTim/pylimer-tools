#include "../../src/pylimer_tools_cpp/utils/DataFileParser.h"
#include "../../src/pylimer_tools_cpp/utils/DataFileParser2.h"
#include "../../src/pylimer_tools_cpp/utils/DumpFileParser.h"
#include "../../src/pylimer_tools_cpp/utils/DumpFileParser2.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include <string>
extern "C"
{
#include <igraph/igraph.h>
}

namespace pu = pylimer_tools::utils;

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

    pu::DumpFileParser2 parser5 =
      pu::DumpFileParser2(suspectedPath + "big_dump_file.lammpstrj");
    CHECK(parser5.getLength() == parser.getLength());

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
}
