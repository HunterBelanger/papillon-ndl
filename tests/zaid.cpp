#include <gtest/gtest.h>

#include <PapillonNDL/element.hpp>
#include <PapillonNDL/isotope.hpp>
#include <PapillonNDL/nuclide.hpp>
#include <PapillonNDL/zaid.hpp>
#include <functional>

namespace pndl {
namespace {

//================================================
// ZAID Tests
TEST(ZAID, ZAID) {
  ZAID zaid(6, 12);

  EXPECT_EQ(zaid.Z(), 6);
  EXPECT_EQ(zaid.A(), 12);
  EXPECT_EQ(zaid.zaid(), 6012);

  ZAID zaid2(6, 13);

  EXPECT_TRUE(zaid == zaid);
  EXPECT_FALSE(zaid == zaid2);

  EXPECT_TRUE(zaid < zaid2);
  EXPECT_FALSE(zaid2 < zaid);

  std::hash<ZAID> zhash;
  std::hash<uint32_t> uhash;
  EXPECT_EQ(zhash(zaid), uhash(zaid.zaid()));
  EXPECT_EQ(zhash(zaid2), uhash(zaid2.zaid()));
}

//================================================
// Element Tests
TEST(Element, Element) {
  Element U(92);
  Element Pu(94);

  ZAID zaid_Th232(90, 232);
  Element Th232(zaid_Th232);

  EXPECT_EQ(U.Z(), 92);
  EXPECT_EQ(Pu.Z(), 94);
  EXPECT_EQ(Th232.Z(), 90);

  EXPECT_EQ(U.atomic_number(), 92);
  EXPECT_EQ(Pu.atomic_number(), 94);
  EXPECT_EQ(Th232.atomic_number(), 90);

  EXPECT_TRUE(U.symbol() == "U");
  EXPECT_TRUE(Pu.symbol() == "Pu");
  EXPECT_TRUE(Th232.symbol() == "Th");

  EXPECT_TRUE(U.name() == "Uranium");
  EXPECT_TRUE(Pu.name() == "Plutonium");
  EXPECT_TRUE(Th232.name() == "Thorium");

  EXPECT_EQ(U.zaid().zaid(), 92000);
  EXPECT_EQ(Pu.zaid().zaid(), 94000);
  EXPECT_EQ(Th232.zaid().zaid(), 90000);

  EXPECT_TRUE(U == U);
  EXPECT_FALSE(U == Pu);

  EXPECT_TRUE(U < Pu);
  EXPECT_FALSE(Pu < U);

  std::hash<Element> ehash;
  std::hash<uint8_t> uhash;
  EXPECT_EQ(ehash(U), uhash(U.Z()));
  EXPECT_EQ(ehash(Pu), uhash(Pu.Z()));

  Element U_from_symbol("U");
  EXPECT_TRUE(U == U_from_symbol);

  Element U_from_name("Uranium");
  EXPECT_TRUE(U == U_from_name);
}

//================================================
// Isotope Tests
TEST(Isotope, Isotope) {
  Element U(92);
  Isotope U235(U, 235);
  EXPECT_EQ(U235.Z(), 92);
  EXPECT_EQ(U235.atomic_number(), 92);
  EXPECT_EQ(U235.A(), 235);
  EXPECT_EQ(U235.atomic_mass(), 235);
  EXPECT_EQ(U235.zaid().zaid(), 92235);
  EXPECT_TRUE(U235.symbol() == "U235");
  EXPECT_TRUE(U235.element_symbol() == "U");
  EXPECT_TRUE(U235.element_name() == "Uranium");

  Isotope Pu239(94, 239);
  EXPECT_EQ(Pu239.Z(), 94);
  EXPECT_EQ(Pu239.atomic_number(), 94);
  EXPECT_EQ(Pu239.A(), 239);
  EXPECT_EQ(Pu239.atomic_mass(), 239);
  EXPECT_EQ(Pu239.zaid().zaid(), 94239);
  EXPECT_TRUE(Pu239.symbol() == "Pu239");
  EXPECT_TRUE(Pu239.element_symbol() == "Pu");
  EXPECT_TRUE(Pu239.element_name() == "Plutonium");

  ZAID zaid_Th232(90, 232);
  Isotope Th232(zaid_Th232);
  EXPECT_EQ(Th232.Z(), 90);
  EXPECT_EQ(Th232.atomic_number(), 90);
  EXPECT_EQ(Th232.A(), 232);
  EXPECT_EQ(Th232.atomic_mass(), 232);
  EXPECT_EQ(Th232.zaid().zaid(), 90232);
  EXPECT_TRUE(Th232.symbol() == "Th232");
  EXPECT_TRUE(Th232.element_symbol() == "Th");
  EXPECT_TRUE(Th232.element_name() == "Thorium");

  EXPECT_TRUE(U235 == U235);
  EXPECT_TRUE(Pu239 == Pu239);
  EXPECT_FALSE(U235 == Pu239);
  Isotope U233(92, 233);
  EXPECT_FALSE(U235 == U233);

  EXPECT_TRUE(U233 < U235);
  EXPECT_FALSE(U235 < U233);

  std::hash<Isotope> ihash;
  std::hash<uint32_t> uhash;
  EXPECT_EQ(ihash(U235), uhash(U235.zaid().zaid()));

  EXPECT_ANY_THROW(Isotope(U, 91));
  EXPECT_NO_THROW(Isotope(U, 92));
  EXPECT_ANY_THROW(Isotope(U, 300));
  EXPECT_NO_THROW(Isotope(U, 299));
  EXPECT_ANY_THROW(Isotope(119, 200));
  EXPECT_ANY_THROW(Isotope(0, 0));
  EXPECT_NO_THROW(Isotope(118, 200));
}

//================================================
// Nuclide Tests
TEST(Nuclide, Nuclide) {
  Element U(92);
  Isotope U235_iso(U, 235);
  Nuclide U235(U235_iso);
  EXPECT_EQ(U235.Z(), 92);
  EXPECT_EQ(U235.atomic_number(), 92);
  EXPECT_EQ(U235.A(), 235);
  EXPECT_EQ(U235.atomic_mass(), 235);
  EXPECT_EQ(U235.zaid().zaid(), 92235);
  EXPECT_TRUE(U235.symbol() == "U235");
  EXPECT_TRUE(U235.isotope_symbol() == "U235");
  EXPECT_TRUE(U235.element_symbol() == "U");
  EXPECT_TRUE(U235.element_name() == "Uranium");

  Nuclide U235m1(92, 235, 1);
  EXPECT_TRUE(U235m1.symbol() == "U235m1");
  EXPECT_TRUE(U235m1.isotope_symbol() == "U235");

  EXPECT_FALSE(U235 == U235m1);
  EXPECT_TRUE(U235 == U235);

  ZAID zaid_U235m2(92, 735);
  Nuclide U235m2(zaid_U235m2);
  EXPECT_EQ(U235m2.Z(), 92);
  EXPECT_EQ(U235m2.atomic_number(), 92);
  EXPECT_EQ(U235m2.A(), 235);
  EXPECT_EQ(U235m2.atomic_mass(), 235);
  EXPECT_EQ(U235m2.zaid().zaid(), 92735);
  EXPECT_TRUE(U235m2.symbol() == "U235m2");
  EXPECT_TRUE(U235m2.isotope_symbol() == "U235");
  EXPECT_TRUE(U235m2.element_symbol() == "U");
  EXPECT_TRUE(U235m2.element_name() == "Uranium");

  Nuclide U236(92, 236);
  Nuclide Pu239(94, 239);
  EXPECT_TRUE(U235 < U235m1);
  EXPECT_FALSE(U235m1 < U235);
  EXPECT_TRUE(U235 < U236);
  EXPECT_TRUE(U235m1 < U236);
  EXPECT_FALSE(U236 < U235m1);
  EXPECT_FALSE(U236 < U235);
  EXPECT_FALSE(Pu239 < U235);
  EXPECT_FALSE(Pu239 < U235m1);
  EXPECT_TRUE(U235 < Pu239);
  EXPECT_TRUE(U235m1 < Pu239);

  Nuclide Pm143("Pm143");
  EXPECT_EQ(Pm143.zaid().zaid(), 61143);

  std::hash<Nuclide> nhash;
  std::hash<uint32_t> uhash;
  EXPECT_EQ(nhash(U235m1), uhash(U235m1.zaid().zaid()));
}

}  // namespace
}  // namespace pndl
