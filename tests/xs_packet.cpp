#include <gtest/gtest.h>

#include <PapillonNDL/xs_packet.hpp>

namespace pndl {
namespace {

//================================================
// XSPacket Tests
TEST(XSPacket, Add) {
  XSPacket xs1;
  xs1.total = 1.;
  xs1.elastic = 2.;
  xs1.inelastic = 3.;
  xs1.absorption = 4.;
  xs1.fission = 5.;
  xs1.capture = 6.;
  xs1.heating = 7.;

  XSPacket xs2 = xs1;

  XSPacket xs3 = xs1 + xs2;

  EXPECT_DOUBLE_EQ(xs3.total, 2.*xs1.total);
  EXPECT_DOUBLE_EQ(xs3.elastic, 2.*xs1.elastic);
  EXPECT_DOUBLE_EQ(xs3.inelastic, 2.*xs1.inelastic);
  EXPECT_DOUBLE_EQ(xs3.absorption, 2.*xs1.absorption);
  EXPECT_DOUBLE_EQ(xs3.fission, 2.*xs1.fission);
  EXPECT_DOUBLE_EQ(xs3.capture, 2.*xs1.capture);
  EXPECT_DOUBLE_EQ(xs3.heating, 2.*xs1.heating);
}

TEST(XSPacket, Sub) {
  XSPacket xs1;
  xs1.total = 1.;
  xs1.elastic = 2.;
  xs1.inelastic = 3.;
  xs1.absorption = 4.;
  xs1.fission = 5.;
  xs1.capture = 6.;
  xs1.heating = 7.;

  XSPacket xs2 = xs1;

  XSPacket xs3 = xs1 - xs2;

  EXPECT_DOUBLE_EQ(xs3.total, 0.);
  EXPECT_DOUBLE_EQ(xs3.elastic, 0.);
  EXPECT_DOUBLE_EQ(xs3.inelastic, 0.);
  EXPECT_DOUBLE_EQ(xs3.absorption, 0.);
  EXPECT_DOUBLE_EQ(xs3.fission, 0.);
  EXPECT_DOUBLE_EQ(xs3.capture, 0.);
  EXPECT_DOUBLE_EQ(xs3.heating, 0.);
}

TEST(XSPacket, Mult) {
  XSPacket xs1;
  xs1.total = 1.;
  xs1.elastic = 2.;
  xs1.inelastic = 3.;
  xs1.absorption = 4.;
  xs1.fission = 5.;
  xs1.capture = 6.;
  xs1.heating = 7.;

  XSPacket xs2 = xs1 * 2.;

  EXPECT_DOUBLE_EQ(xs2.total, 2.*xs1.total);
  EXPECT_DOUBLE_EQ(xs2.elastic, 2.*xs1.elastic);
  EXPECT_DOUBLE_EQ(xs2.inelastic, 2.*xs1.inelastic);
  EXPECT_DOUBLE_EQ(xs2.absorption, 2.*xs1.absorption);
  EXPECT_DOUBLE_EQ(xs2.fission, 2.*xs1.fission);
  EXPECT_DOUBLE_EQ(xs2.capture, 2.*xs1.capture);
  EXPECT_DOUBLE_EQ(xs2.heating, 2.*xs1.heating);

  XSPacket xs3 = 3. * xs1;

  EXPECT_DOUBLE_EQ(xs3.total, 3.*xs1.total);
  EXPECT_DOUBLE_EQ(xs3.elastic, 3.*xs1.elastic);
  EXPECT_DOUBLE_EQ(xs3.inelastic, 3.*xs1.inelastic);
  EXPECT_DOUBLE_EQ(xs3.absorption, 3.*xs1.absorption);
  EXPECT_DOUBLE_EQ(xs3.fission, 3.*xs1.fission);
  EXPECT_DOUBLE_EQ(xs3.capture, 3.*xs1.capture);
  EXPECT_DOUBLE_EQ(xs3.heating, 3.*xs1.heating);
}

TEST(XSPacket, Div) {
  XSPacket xs1;
  xs1.total = 1.;
  xs1.elastic = 2.;
  xs1.inelastic = 3.;
  xs1.absorption = 4.;
  xs1.fission = 5.;
  xs1.capture = 6.;
  xs1.heating = 7.;

  XSPacket xs2 = xs1 / 2.;

  EXPECT_DOUBLE_EQ(xs2.total, 0.5*xs1.total);
  EXPECT_DOUBLE_EQ(xs2.elastic, 0.5*xs1.elastic);
  EXPECT_DOUBLE_EQ(xs2.inelastic, 0.5*xs1.inelastic);
  EXPECT_DOUBLE_EQ(xs2.absorption, 0.5*xs1.absorption);
  EXPECT_DOUBLE_EQ(xs2.fission, 0.5*xs1.fission);
  EXPECT_DOUBLE_EQ(xs2.capture, 0.5*xs1.capture);
  EXPECT_DOUBLE_EQ(xs2.heating, 0.5*xs1.heating);
}

TEST(XSPacket, AddAssign) {
  XSPacket xs1;
  xs1.total = 1.;
  xs1.elastic = 2.;
  xs1.inelastic = 3.;
  xs1.absorption = 4.;
  xs1.fission = 5.;
  xs1.capture = 6.;
  xs1.heating = 7.;

  XSPacket xs2 = xs1;

  xs1 += xs2;

  EXPECT_DOUBLE_EQ(xs1.total, 2.*xs2.total);
  EXPECT_DOUBLE_EQ(xs1.elastic, 2.*xs2.elastic);
  EXPECT_DOUBLE_EQ(xs1.inelastic, 2.*xs2.inelastic);
  EXPECT_DOUBLE_EQ(xs1.absorption, 2.*xs2.absorption);
  EXPECT_DOUBLE_EQ(xs1.fission, 2.*xs2.fission);
  EXPECT_DOUBLE_EQ(xs1.capture, 2.*xs2.capture);
  EXPECT_DOUBLE_EQ(xs1.heating, 2.*xs2.heating);
}

TEST(XSPacket, SubAssign) {
  XSPacket xs1;
  xs1.total = 1.;
  xs1.elastic = 2.;
  xs1.inelastic = 3.;
  xs1.absorption = 4.;
  xs1.fission = 5.;
  xs1.capture = 6.;
  xs1.heating = 7.;

  XSPacket xs2 = xs1;

  xs1 -= xs2;

  EXPECT_DOUBLE_EQ(xs1.total, 0.);
  EXPECT_DOUBLE_EQ(xs1.elastic, 0.);
  EXPECT_DOUBLE_EQ(xs1.inelastic, 0.);
  EXPECT_DOUBLE_EQ(xs1.absorption, 0.);
  EXPECT_DOUBLE_EQ(xs1.fission, 0.);
  EXPECT_DOUBLE_EQ(xs1.capture, 0.);
  EXPECT_DOUBLE_EQ(xs1.heating, 0.);
}

TEST(XSPacket, MultAssign) {
  XSPacket xs1;
  xs1.total = 1.;
  xs1.elastic = 2.;
  xs1.inelastic = 3.;
  xs1.absorption = 4.;
  xs1.fission = 5.;
  xs1.capture = 6.;
  xs1.heating = 7.;

  XSPacket xs2 = xs1;
  xs2 *= 2.;

  EXPECT_DOUBLE_EQ(xs2.total, 2.*xs1.total);
  EXPECT_DOUBLE_EQ(xs2.elastic, 2.*xs1.elastic);
  EXPECT_DOUBLE_EQ(xs2.inelastic, 2.*xs1.inelastic);
  EXPECT_DOUBLE_EQ(xs2.absorption, 2.*xs1.absorption);
  EXPECT_DOUBLE_EQ(xs2.fission, 2.*xs1.fission);
  EXPECT_DOUBLE_EQ(xs2.capture, 2.*xs1.capture);
  EXPECT_DOUBLE_EQ(xs2.heating, 2.*xs1.heating);
}

TEST(XSPacket, DivAssign) {
  XSPacket xs1;
  xs1.total = 1.;
  xs1.elastic = 2.;
  xs1.inelastic = 3.;
  xs1.absorption = 4.;
  xs1.fission = 5.;
  xs1.capture = 6.;
  xs1.heating = 7.;

  XSPacket xs2 = xs1;
  xs2 /= 2.;

  EXPECT_DOUBLE_EQ(xs2.total, 0.5*xs1.total);
  EXPECT_DOUBLE_EQ(xs2.elastic, 0.5*xs1.elastic);
  EXPECT_DOUBLE_EQ(xs2.inelastic, 0.5*xs1.inelastic);
  EXPECT_DOUBLE_EQ(xs2.absorption, 0.5*xs1.absorption);
  EXPECT_DOUBLE_EQ(xs2.fission, 0.5*xs1.fission);
  EXPECT_DOUBLE_EQ(xs2.capture, 0.5*xs1.capture);
  EXPECT_DOUBLE_EQ(xs2.heating, 0.5*xs1.heating);
}

}
}

