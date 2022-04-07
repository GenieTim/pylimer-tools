#include "../../src/pylimer_tools_cpp/utils/DataFileParser.h"
#include "../../src/pylimer_tools_cpp/utils/DumpFileParser.h"
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include <string>
extern "C"
{
#include <igraph/igraph.h>
}

namespace pu = pylimer_tools::utils;

TEST_CASE("DumpFileParser can be used", "[utils][DumpFileParser]")
{
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));

  SECTION("Reading from dump file works")
  {
    pu::DumpFileParser parser =
      pu::DumpFileParser(suspectedPath + "lammps_dump_small.lammpstrj");
    REQUIRE(parser.getLength() == 1);
    REQUIRE_NOTHROW(parser.read());
    REQUIRE(parser.getLength() == 1);
    REQUIRE_THROWS(parser.readNGroups(1, 10));
    REQUIRE_THROWS(parser.getValuesForAt<int>(0, "BOX BOUNDS", "notexisting"));
    // call copy constructor
    pu::DumpFileParser parser2 = parser;
    REQUIRE(parser2.getLength() == 1);

    // test other reading capabilities
    pu::DumpFileParser parser3 =
      pu::DumpFileParser(suspectedPath + "lammps_dump_small.lammpstrj");
    REQUIRE(parser.getValuesForAt<double>(0, "BOX BOUNDS", 1).size() == 3);
    REQUIRE_THROWS(parser.getValuesForAt<double>(0, "NOT EXISTING", 9));
  }

  SECTION("Reading from data files works")
  {
    pu::DataFileParser parser = pu::DataFileParser();
    parser.read(suspectedPath + "lammps_data_file.out");
    REQUIRE(parser.getNrOfAtoms() == 3000);
    // call copy constructor
    pu::DataFileParser parser2 = parser;
    REQUIRE(parser2.getNrOfAtoms() == 3000);
  }
}
