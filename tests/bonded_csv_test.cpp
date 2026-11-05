//
// Copyright 2026 Gabriel Lando
//
// Comma-separated block parameters (serials and per-radio values).
//

#include "bonded_csv.hpp"
#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <string>
#include <vector>

using namespace bonded::detail;

BOOST_AUTO_TEST_CASE(test_split_csv_trims_and_drops_empty)
{
    const std::vector<std::string> expect = {"SERIAL_A", "SERIAL_B"};
    BOOST_CHECK(split_csv("SERIAL_A,SERIAL_B") == expect);
    BOOST_CHECK(split_csv("  SERIAL_A , \tSERIAL_B  ") == expect);
    BOOST_CHECK(split_csv(",SERIAL_A,,SERIAL_B,") == expect);
    BOOST_CHECK(split_csv("").empty());
    BOOST_CHECK(split_csv(" , ").empty());
    BOOST_CHECK(split_csv("one") == std::vector<std::string>{"one"});
}

BOOST_AUTO_TEST_CASE(test_parse_csv_doubles)
{
    const auto v = parse_csv_doubles("2.40625e9, 2.41375e9", "freq_plan_csv");
    BOOST_REQUIRE_EQUAL(v.size(), size_t(2));
    BOOST_CHECK_EQUAL(v[0], 2.40625e9);
    BOOST_CHECK_EQUAL(v[1], 2.41375e9);
    BOOST_CHECK(parse_csv_doubles("", "x").empty());
    BOOST_CHECK_EQUAL(parse_csv_doubles("-0.5", "x")[0], -0.5);
}

BOOST_AUTO_TEST_CASE(test_parse_csv_doubles_rejects_garbage)
{
    BOOST_CHECK_THROW(parse_csv_doubles("1e6,abc", "freq_plan_csv"), std::invalid_argument);
    BOOST_CHECK_THROW(parse_csv_doubles("1e6x", "freq_plan_csv"), std::invalid_argument);
    BOOST_CHECK_THROW(parse_csv_doubles("1 2", "freq_plan_csv"), std::invalid_argument);
    try {
        parse_csv_doubles("0.1,zz", "gain_plan_csv");
        BOOST_FAIL("expected an exception");
    } catch (const std::invalid_argument& e) {
        BOOST_CHECK(std::string(e.what()).find("gain_plan_csv") != std::string::npos);
        BOOST_CHECK(std::string(e.what()).find("zz") != std::string::npos);
    }
}
