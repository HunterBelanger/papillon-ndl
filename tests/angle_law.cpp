#include <PapillonNDL/isotropic.hpp>
#include <PapillonNDL/equiprobable_angle_bins.hpp>
#include <PapillonNDL/angle_table.hpp>
#include <vector>
#include <gtest/gtest.h>

namespace pndl {
namespace {

//==============================================================================
// Isotropic Tests
TEST(Isotropic, SampleMu) {
  Isotropic iso;

  const std::vector<double> xi {0., 0.25, 0.5, 0.75, 1.};
  std::size_t xi_i = 0;
  auto rng = [&xi, &xi_i]() {
    return xi[xi_i++]; 
  };

  EXPECT_DOUBLE_EQ(iso.sample_mu(rng), -1.);
  EXPECT_DOUBLE_EQ(iso.sample_mu(rng), -0.5);
  EXPECT_DOUBLE_EQ(iso.sample_mu(rng), 0.);
  EXPECT_DOUBLE_EQ(iso.sample_mu(rng), 0.5);
  EXPECT_DOUBLE_EQ(iso.sample_mu(rng), 1.);
}

TEST(Isotropic, PDF) {
  Isotropic iso;

  EXPECT_DOUBLE_EQ(iso.pdf(-1.), 0.5);
  EXPECT_DOUBLE_EQ(iso.pdf(-0.5), 0.5);
  EXPECT_DOUBLE_EQ(iso.pdf(0.), 0.5);
  EXPECT_DOUBLE_EQ(iso.pdf(0.5), 0.5);
  EXPECT_DOUBLE_EQ(iso.pdf(1.), 0.5);
}

//==============================================================================
// EquiprobableAngleBins Tests
TEST(EquiprobableAngleBins, SampleMu) {
  // These bounds correspoond with an isotropic distribution.
  std::vector<double> bounds
      {-1.    , -0.9375, -0.875 , -0.8125, -0.75  , -0.6875, -0.625 ,
       -0.5625, -0.5   , -0.4375, -0.375 , -0.3125, -0.25  , -0.1875,
       -0.125 , -0.0625,  0.    ,  0.0625,  0.125 ,  0.1875,  0.25  ,
        0.3125,  0.375 ,  0.4375,  0.5   ,  0.5625,  0.625 ,  0.6875,
        0.75  ,  0.8125,  0.875 ,  0.9375,  1.};
  EquiprobableAngleBins bins(bounds);

  const std::vector<double> xi {0., 0.25, 0.5, 0.75, 1.};
  std::size_t xi_i = 0;
  auto rng = [&xi, &xi_i]() {
    return xi[xi_i++]; 
  };

  EXPECT_DOUBLE_EQ(bins.sample_mu(rng), -1.);
  EXPECT_DOUBLE_EQ(bins.sample_mu(rng), -0.5);
  EXPECT_DOUBLE_EQ(bins.sample_mu(rng), 0.);
  EXPECT_DOUBLE_EQ(bins.sample_mu(rng), 0.5);
  EXPECT_DOUBLE_EQ(bins.sample_mu(rng), 1.);
}

TEST(EquiprobableAngleBins, PDF) {
  // These bounds correspoond with an isotropic distribution.
  std::vector<double> bounds
      {-1.    , -0.9375, -0.875 , -0.8125, -0.75  , -0.6875, -0.625 ,
       -0.5625, -0.5   , -0.4375, -0.375 , -0.3125, -0.25  , -0.1875,
       -0.125 , -0.0625,  0.    ,  0.0625,  0.125 ,  0.1875,  0.25  ,
        0.3125,  0.375 ,  0.4375,  0.5   ,  0.5625,  0.625 ,  0.6875,
        0.75  ,  0.8125,  0.875 ,  0.9375,  1.};
  EquiprobableAngleBins bins(bounds);
  
  EXPECT_DOUBLE_EQ(bins.pdf(-1.), 0.5);
  EXPECT_DOUBLE_EQ(bins.pdf(-0.5), 0.5);
  EXPECT_DOUBLE_EQ(bins.pdf(0.), 0.5);
  EXPECT_DOUBLE_EQ(bins.pdf(0.5), 0.5);
  EXPECT_DOUBLE_EQ(bins.pdf(1.), 0.5);
}

TEST(EquiprobableAngleBins, Size) {
  // These bounds correspoond with an isotropic distribution.
  std::vector<double> bounds
      {-1.    , -0.9375, -0.875 , -0.8125, -0.75  , -0.6875, -0.625 ,
       -0.5625, -0.5   , -0.4375, -0.375 , -0.3125, -0.25  , -0.1875,
       -0.125 , -0.0625,  0.    ,  0.0625,  0.125 ,  0.1875,  0.25  ,
        0.3125,  0.375 ,  0.4375,  0.5   ,  0.5625,  0.625 ,  0.6875,
        0.75  ,  0.8125,  0.875 ,  0.9375,  1.};
  EquiprobableAngleBins bins(bounds);
  
  EXPECT_EQ(bins.size(), 33);
}

TEST(EquiprobableAngleBins, Bounds) {
  // These bounds correspoond with an isotropic distribution.
  std::vector<double> bounds
      {-1.    , -0.9375, -0.875 , -0.8125, -0.75  , -0.6875, -0.625 ,
       -0.5625, -0.5   , -0.4375, -0.375 , -0.3125, -0.25  , -0.1875,
       -0.125 , -0.0625,  0.    ,  0.0625,  0.125 ,  0.1875,  0.25  ,
        0.3125,  0.375 ,  0.4375,  0.5   ,  0.5625,  0.625 ,  0.6875,
        0.75  ,  0.8125,  0.875 ,  0.9375,  1.};
  EquiprobableAngleBins bins(bounds);

  const auto& bin_bounds = bins.bin_bounds();

  ASSERT_EQ(bin_bounds.size(), bounds.size());

  for (std::size_t i = 0; i < bin_bounds.size(); i++) {
    EXPECT_DOUBLE_EQ(bin_bounds[i], bounds[i]); 
  }
}

//==============================================================================
// AngleTable Tests
TEST(AngleTable, Construction) {
  std::vector<double> vals1 {-1.1, 0., 1.};
  std::vector<double> vals2 {-1., 0., 1.1};
  std::vector<double> vals3 {-1., 0., 1.};
  std::vector<double> pdf   {0.5, 0.5, 0.5};
  std::vector<double> cdf   {0., 0.5, 1.};

  EXPECT_THROW(AngleTable(vals1,pdf,cdf,Interpolation::LinLin), PNDLException);
  EXPECT_THROW(AngleTable(vals2,pdf,cdf,Interpolation::LinLin), PNDLException);
  EXPECT_NO_THROW(AngleTable(vals3,pdf,cdf,Interpolation::LinLin));
}

TEST(AngleTable, SampleMu) {
  std::vector<double> vals {-1., 0., 1.};
  std::vector<double> pdf   {0.5, 0.5, 0.5};
  std::vector<double> cdf   {0., 0.5, 1.};
  AngleTable tab(vals, pdf, cdf, Interpolation::LinLin);
  
  const std::vector<double> xi {0., 0.25, 0.5, 0.75, 1.};
  std::size_t xi_i = 0;
  auto rng = [&xi, &xi_i]() {
    return xi[xi_i++]; 
  };

  EXPECT_DOUBLE_EQ(tab.sample_mu(rng), -1.);
  EXPECT_DOUBLE_EQ(tab.sample_mu(rng), -0.5);
  EXPECT_DOUBLE_EQ(tab.sample_mu(rng), 0.);
  EXPECT_DOUBLE_EQ(tab.sample_mu(rng), 0.5);
  EXPECT_DOUBLE_EQ(tab.sample_mu(rng), 1.);
}

TEST(AngleTable, PDF) {
  std::vector<double> vals {-1., 0., 1.};
  std::vector<double> pdf   {0.5, 0.5, 0.5};
  std::vector<double> cdf   {0., 0.5, 1.};
  AngleTable tab(vals, pdf, cdf, Interpolation::LinLin);

  EXPECT_DOUBLE_EQ(tab.pdf(-1.), 0.5);
  EXPECT_DOUBLE_EQ(tab.pdf(-0.5), 0.5);
  EXPECT_DOUBLE_EQ(tab.pdf(0.), 0.5);
  EXPECT_DOUBLE_EQ(tab.pdf(0.5), 0.5);
  EXPECT_DOUBLE_EQ(tab.pdf(1.), 0.5);
}

TEST(AngleTable, Size) {
  std::vector<double> cosines {-1., -0.25, 0.25, 1.};
  std::vector<double> pdf {0., 0.25, 0.25, 1.};
  std::vector<double> cdf {0., 0.125, 0.375, 1.};
  AngleTable tab(cosines, pdf, cdf, Interpolation::LinLin);

  EXPECT_EQ(tab.size(), cosines.size());
}

TEST(AngleTable, ValuesGrid) {
  std::vector<double> cosines {-1., -0.25, 0.25, 1.};
  std::vector<double> pdf {0., 0.25, 0.25, 1.};
  std::vector<double> cdf {0., 0.125, 0.375, 1.};
  AngleTable tab(cosines, pdf, cdf, Interpolation::LinLin);

  const auto& vals = tab.cosines();

  ASSERT_EQ(vals.size(), tab.size());

  for (std::size_t i = 0; i < vals.size(); i++) {
    EXPECT_DOUBLE_EQ(vals[i], cosines[i]);
  }
}

TEST(AngleTable, PDFGrid) {
  std::vector<double> cosines {-1., -0.25, 0.25, 1.};
  std::vector<double> pdf {0., 0.25, 0.25, 1.};
  std::vector<double> cdf {0., 0.125, 0.375, 1.};
  AngleTable tab(cosines, pdf, cdf, Interpolation::LinLin);

  const auto& vals = tab.pdf();

  ASSERT_EQ(vals.size(), pdf.size());

  for (std::size_t i = 0; i < vals.size(); i++) {
    EXPECT_DOUBLE_EQ(vals[i], pdf[i]);
  }
}

TEST(AngleTable, CDFGrid) {
  std::vector<double> cosines {-1., -0.25, 0.25, 1.};
  std::vector<double> pdf {0., 0.25, 0.25, 1.};
  std::vector<double> cdf {0., 0.125, 0.375, 1.};
  AngleTable tab(cosines, pdf, cdf, Interpolation::LinLin);

  const auto& vals = tab.cdf();

  ASSERT_EQ(vals.size(), cdf.size());

  for (std::size_t i = 0; i < vals.size(); i++) {
    EXPECT_DOUBLE_EQ(vals[i], cdf[i]);
  }
}

TEST(AngleTable, Interpolation) {
  std::vector<double> cosines {-1., -0.25, 0.25, 1.};
  std::vector<double> pdf {0., 0.25, 0.25, 1.};
  std::vector<double> cdf {0., 0.125, 0.375, 1.};
  AngleTable lin(cosines, pdf, cdf, Interpolation::LinLin);
  
  EXPECT_EQ(Interpolation::LinLin, lin.interpolation());

  AngleTable hist(cosines, pdf, cdf, Interpolation::Histogram);
  EXPECT_EQ(Interpolation::Histogram, hist.interpolation());
}

}
}
