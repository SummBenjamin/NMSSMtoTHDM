// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Mon 20 Mar 2017 08:52:43

#include "NMSSMtower_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

NMSSMtower_physical::NMSSMtower_physical()
   :
    MVG(0), MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MSd(Eigen::Array<
       double,6,1>::Zero()), MSv(Eigen::Array<double,3,1>::Zero()), MSu(
       Eigen::Array<double,6,1>::Zero()), MSe(Eigen::Array<double,6,1>::Zero()),
       Mhh(Eigen::Array<double,3,1>::Zero()), MAh(Eigen::Array<double,3,1>::Zero()
       ), MHpm(Eigen::Array<double,2,1>::Zero()), MChi(Eigen::Array<double,5,1>
       ::Zero()), MCha(Eigen::Array<double,2,1>::Zero()), MFe(Eigen::Array<double,
       3,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<
       double,3,1>::Zero()), MVWm(0), MVP(0), MVZ(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,3,3>::Zero()), ZA(Eigen::Matrix<double,3,
      3>::Zero()), ZP(Eigen::Matrix<double,2,2>::Zero()), ZN(Eigen::Matrix<
      std::complex<double>,5,5>::Zero()), UM(Eigen::Matrix<std::complex<double>,2,
      2>::Zero()), UP(Eigen::Matrix<std::complex<double>,2,2>::Zero()), ZEL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZER(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZDL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZDR(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZZ(Eigen::Matrix<double,2,2>::Zero())

{
}

void NMSSMtower_physical::clear()
{
   MVG = 0.;
   MGlu = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,5,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,5,5>::Zero();
   MCha = Eigen::Matrix<double,2,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   ZEL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZER = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   ZDL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   ZUL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZUR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;

}

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void NMSSMtower_physical::convert_to_hk()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void NMSSMtower_physical::convert_to_slha()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

Eigen::ArrayXd NMSSMtower_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

   pars.conservativeResize(370);

   pars(53) = ZD(0,0);
   pars(54) = ZD(0,1);
   pars(55) = ZD(0,2);
   pars(56) = ZD(0,3);
   pars(57) = ZD(0,4);
   pars(58) = ZD(0,5);
   pars(59) = ZD(1,0);
   pars(60) = ZD(1,1);
   pars(61) = ZD(1,2);
   pars(62) = ZD(1,3);
   pars(63) = ZD(1,4);
   pars(64) = ZD(1,5);
   pars(65) = ZD(2,0);
   pars(66) = ZD(2,1);
   pars(67) = ZD(2,2);
   pars(68) = ZD(2,3);
   pars(69) = ZD(2,4);
   pars(70) = ZD(2,5);
   pars(71) = ZD(3,0);
   pars(72) = ZD(3,1);
   pars(73) = ZD(3,2);
   pars(74) = ZD(3,3);
   pars(75) = ZD(3,4);
   pars(76) = ZD(3,5);
   pars(77) = ZD(4,0);
   pars(78) = ZD(4,1);
   pars(79) = ZD(4,2);
   pars(80) = ZD(4,3);
   pars(81) = ZD(4,4);
   pars(82) = ZD(4,5);
   pars(83) = ZD(5,0);
   pars(84) = ZD(5,1);
   pars(85) = ZD(5,2);
   pars(86) = ZD(5,3);
   pars(87) = ZD(5,4);
   pars(88) = ZD(5,5);
   pars(89) = ZV(0,0);
   pars(90) = ZV(0,1);
   pars(91) = ZV(0,2);
   pars(92) = ZV(1,0);
   pars(93) = ZV(1,1);
   pars(94) = ZV(1,2);
   pars(95) = ZV(2,0);
   pars(96) = ZV(2,1);
   pars(97) = ZV(2,2);
   pars(98) = ZU(0,0);
   pars(99) = ZU(0,1);
   pars(100) = ZU(0,2);
   pars(101) = ZU(0,3);
   pars(102) = ZU(0,4);
   pars(103) = ZU(0,5);
   pars(104) = ZU(1,0);
   pars(105) = ZU(1,1);
   pars(106) = ZU(1,2);
   pars(107) = ZU(1,3);
   pars(108) = ZU(1,4);
   pars(109) = ZU(1,5);
   pars(110) = ZU(2,0);
   pars(111) = ZU(2,1);
   pars(112) = ZU(2,2);
   pars(113) = ZU(2,3);
   pars(114) = ZU(2,4);
   pars(115) = ZU(2,5);
   pars(116) = ZU(3,0);
   pars(117) = ZU(3,1);
   pars(118) = ZU(3,2);
   pars(119) = ZU(3,3);
   pars(120) = ZU(3,4);
   pars(121) = ZU(3,5);
   pars(122) = ZU(4,0);
   pars(123) = ZU(4,1);
   pars(124) = ZU(4,2);
   pars(125) = ZU(4,3);
   pars(126) = ZU(4,4);
   pars(127) = ZU(4,5);
   pars(128) = ZU(5,0);
   pars(129) = ZU(5,1);
   pars(130) = ZU(5,2);
   pars(131) = ZU(5,3);
   pars(132) = ZU(5,4);
   pars(133) = ZU(5,5);
   pars(134) = ZE(0,0);
   pars(135) = ZE(0,1);
   pars(136) = ZE(0,2);
   pars(137) = ZE(0,3);
   pars(138) = ZE(0,4);
   pars(139) = ZE(0,5);
   pars(140) = ZE(1,0);
   pars(141) = ZE(1,1);
   pars(142) = ZE(1,2);
   pars(143) = ZE(1,3);
   pars(144) = ZE(1,4);
   pars(145) = ZE(1,5);
   pars(146) = ZE(2,0);
   pars(147) = ZE(2,1);
   pars(148) = ZE(2,2);
   pars(149) = ZE(2,3);
   pars(150) = ZE(2,4);
   pars(151) = ZE(2,5);
   pars(152) = ZE(3,0);
   pars(153) = ZE(3,1);
   pars(154) = ZE(3,2);
   pars(155) = ZE(3,3);
   pars(156) = ZE(3,4);
   pars(157) = ZE(3,5);
   pars(158) = ZE(4,0);
   pars(159) = ZE(4,1);
   pars(160) = ZE(4,2);
   pars(161) = ZE(4,3);
   pars(162) = ZE(4,4);
   pars(163) = ZE(4,5);
   pars(164) = ZE(5,0);
   pars(165) = ZE(5,1);
   pars(166) = ZE(5,2);
   pars(167) = ZE(5,3);
   pars(168) = ZE(5,4);
   pars(169) = ZE(5,5);
   pars(170) = ZH(0,0);
   pars(171) = ZH(0,1);
   pars(172) = ZH(0,2);
   pars(173) = ZH(1,0);
   pars(174) = ZH(1,1);
   pars(175) = ZH(1,2);
   pars(176) = ZH(2,0);
   pars(177) = ZH(2,1);
   pars(178) = ZH(2,2);
   pars(179) = ZA(0,0);
   pars(180) = ZA(0,1);
   pars(181) = ZA(0,2);
   pars(182) = ZA(1,0);
   pars(183) = ZA(1,1);
   pars(184) = ZA(1,2);
   pars(185) = ZA(2,0);
   pars(186) = ZA(2,1);
   pars(187) = ZA(2,2);
   pars(188) = ZP(0,0);
   pars(189) = ZP(0,1);
   pars(190) = ZP(1,0);
   pars(191) = ZP(1,1);
   pars(192) = Re(ZN(0,0));
   pars(193) = Im(ZN(0,0));
   pars(194) = Re(ZN(0,1));
   pars(195) = Im(ZN(0,1));
   pars(196) = Re(ZN(0,2));
   pars(197) = Im(ZN(0,2));
   pars(198) = Re(ZN(0,3));
   pars(199) = Im(ZN(0,3));
   pars(200) = Re(ZN(0,4));
   pars(201) = Im(ZN(0,4));
   pars(202) = Re(ZN(1,0));
   pars(203) = Im(ZN(1,0));
   pars(204) = Re(ZN(1,1));
   pars(205) = Im(ZN(1,1));
   pars(206) = Re(ZN(1,2));
   pars(207) = Im(ZN(1,2));
   pars(208) = Re(ZN(1,3));
   pars(209) = Im(ZN(1,3));
   pars(210) = Re(ZN(1,4));
   pars(211) = Im(ZN(1,4));
   pars(212) = Re(ZN(2,0));
   pars(213) = Im(ZN(2,0));
   pars(214) = Re(ZN(2,1));
   pars(215) = Im(ZN(2,1));
   pars(216) = Re(ZN(2,2));
   pars(217) = Im(ZN(2,2));
   pars(218) = Re(ZN(2,3));
   pars(219) = Im(ZN(2,3));
   pars(220) = Re(ZN(2,4));
   pars(221) = Im(ZN(2,4));
   pars(222) = Re(ZN(3,0));
   pars(223) = Im(ZN(3,0));
   pars(224) = Re(ZN(3,1));
   pars(225) = Im(ZN(3,1));
   pars(226) = Re(ZN(3,2));
   pars(227) = Im(ZN(3,2));
   pars(228) = Re(ZN(3,3));
   pars(229) = Im(ZN(3,3));
   pars(230) = Re(ZN(3,4));
   pars(231) = Im(ZN(3,4));
   pars(232) = Re(ZN(4,0));
   pars(233) = Im(ZN(4,0));
   pars(234) = Re(ZN(4,1));
   pars(235) = Im(ZN(4,1));
   pars(236) = Re(ZN(4,2));
   pars(237) = Im(ZN(4,2));
   pars(238) = Re(ZN(4,3));
   pars(239) = Im(ZN(4,3));
   pars(240) = Re(ZN(4,4));
   pars(241) = Im(ZN(4,4));
   pars(242) = Re(UM(0,0));
   pars(243) = Im(UM(0,0));
   pars(244) = Re(UM(0,1));
   pars(245) = Im(UM(0,1));
   pars(246) = Re(UM(1,0));
   pars(247) = Im(UM(1,0));
   pars(248) = Re(UM(1,1));
   pars(249) = Im(UM(1,1));
   pars(250) = Re(UP(0,0));
   pars(251) = Im(UP(0,0));
   pars(252) = Re(UP(0,1));
   pars(253) = Im(UP(0,1));
   pars(254) = Re(UP(1,0));
   pars(255) = Im(UP(1,0));
   pars(256) = Re(UP(1,1));
   pars(257) = Im(UP(1,1));
   pars(258) = Re(ZEL(0,0));
   pars(259) = Im(ZEL(0,0));
   pars(260) = Re(ZEL(0,1));
   pars(261) = Im(ZEL(0,1));
   pars(262) = Re(ZEL(0,2));
   pars(263) = Im(ZEL(0,2));
   pars(264) = Re(ZEL(1,0));
   pars(265) = Im(ZEL(1,0));
   pars(266) = Re(ZEL(1,1));
   pars(267) = Im(ZEL(1,1));
   pars(268) = Re(ZEL(1,2));
   pars(269) = Im(ZEL(1,2));
   pars(270) = Re(ZEL(2,0));
   pars(271) = Im(ZEL(2,0));
   pars(272) = Re(ZEL(2,1));
   pars(273) = Im(ZEL(2,1));
   pars(274) = Re(ZEL(2,2));
   pars(275) = Im(ZEL(2,2));
   pars(276) = Re(ZER(0,0));
   pars(277) = Im(ZER(0,0));
   pars(278) = Re(ZER(0,1));
   pars(279) = Im(ZER(0,1));
   pars(280) = Re(ZER(0,2));
   pars(281) = Im(ZER(0,2));
   pars(282) = Re(ZER(1,0));
   pars(283) = Im(ZER(1,0));
   pars(284) = Re(ZER(1,1));
   pars(285) = Im(ZER(1,1));
   pars(286) = Re(ZER(1,2));
   pars(287) = Im(ZER(1,2));
   pars(288) = Re(ZER(2,0));
   pars(289) = Im(ZER(2,0));
   pars(290) = Re(ZER(2,1));
   pars(291) = Im(ZER(2,1));
   pars(292) = Re(ZER(2,2));
   pars(293) = Im(ZER(2,2));
   pars(294) = Re(ZDL(0,0));
   pars(295) = Im(ZDL(0,0));
   pars(296) = Re(ZDL(0,1));
   pars(297) = Im(ZDL(0,1));
   pars(298) = Re(ZDL(0,2));
   pars(299) = Im(ZDL(0,2));
   pars(300) = Re(ZDL(1,0));
   pars(301) = Im(ZDL(1,0));
   pars(302) = Re(ZDL(1,1));
   pars(303) = Im(ZDL(1,1));
   pars(304) = Re(ZDL(1,2));
   pars(305) = Im(ZDL(1,2));
   pars(306) = Re(ZDL(2,0));
   pars(307) = Im(ZDL(2,0));
   pars(308) = Re(ZDL(2,1));
   pars(309) = Im(ZDL(2,1));
   pars(310) = Re(ZDL(2,2));
   pars(311) = Im(ZDL(2,2));
   pars(312) = Re(ZDR(0,0));
   pars(313) = Im(ZDR(0,0));
   pars(314) = Re(ZDR(0,1));
   pars(315) = Im(ZDR(0,1));
   pars(316) = Re(ZDR(0,2));
   pars(317) = Im(ZDR(0,2));
   pars(318) = Re(ZDR(1,0));
   pars(319) = Im(ZDR(1,0));
   pars(320) = Re(ZDR(1,1));
   pars(321) = Im(ZDR(1,1));
   pars(322) = Re(ZDR(1,2));
   pars(323) = Im(ZDR(1,2));
   pars(324) = Re(ZDR(2,0));
   pars(325) = Im(ZDR(2,0));
   pars(326) = Re(ZDR(2,1));
   pars(327) = Im(ZDR(2,1));
   pars(328) = Re(ZDR(2,2));
   pars(329) = Im(ZDR(2,2));
   pars(330) = Re(ZUL(0,0));
   pars(331) = Im(ZUL(0,0));
   pars(332) = Re(ZUL(0,1));
   pars(333) = Im(ZUL(0,1));
   pars(334) = Re(ZUL(0,2));
   pars(335) = Im(ZUL(0,2));
   pars(336) = Re(ZUL(1,0));
   pars(337) = Im(ZUL(1,0));
   pars(338) = Re(ZUL(1,1));
   pars(339) = Im(ZUL(1,1));
   pars(340) = Re(ZUL(1,2));
   pars(341) = Im(ZUL(1,2));
   pars(342) = Re(ZUL(2,0));
   pars(343) = Im(ZUL(2,0));
   pars(344) = Re(ZUL(2,1));
   pars(345) = Im(ZUL(2,1));
   pars(346) = Re(ZUL(2,2));
   pars(347) = Im(ZUL(2,2));
   pars(348) = Re(ZUR(0,0));
   pars(349) = Im(ZUR(0,0));
   pars(350) = Re(ZUR(0,1));
   pars(351) = Im(ZUR(0,1));
   pars(352) = Re(ZUR(0,2));
   pars(353) = Im(ZUR(0,2));
   pars(354) = Re(ZUR(1,0));
   pars(355) = Im(ZUR(1,0));
   pars(356) = Re(ZUR(1,1));
   pars(357) = Im(ZUR(1,1));
   pars(358) = Re(ZUR(1,2));
   pars(359) = Im(ZUR(1,2));
   pars(360) = Re(ZUR(2,0));
   pars(361) = Im(ZUR(2,0));
   pars(362) = Re(ZUR(2,1));
   pars(363) = Im(ZUR(2,1));
   pars(364) = Re(ZUR(2,2));
   pars(365) = Im(ZUR(2,2));
   pars(366) = ZZ(0,0);
   pars(367) = ZZ(0,1);
   pars(368) = ZZ(1,0);
   pars(369) = ZZ(1,1);


   return pars;
}

void NMSSMtower_physical::set(const Eigen::ArrayXd& pars)
{
   set_masses(pars);

   ZD(0,0) = pars(53);
   ZD(0,1) = pars(54);
   ZD(0,2) = pars(55);
   ZD(0,3) = pars(56);
   ZD(0,4) = pars(57);
   ZD(0,5) = pars(58);
   ZD(1,0) = pars(59);
   ZD(1,1) = pars(60);
   ZD(1,2) = pars(61);
   ZD(1,3) = pars(62);
   ZD(1,4) = pars(63);
   ZD(1,5) = pars(64);
   ZD(2,0) = pars(65);
   ZD(2,1) = pars(66);
   ZD(2,2) = pars(67);
   ZD(2,3) = pars(68);
   ZD(2,4) = pars(69);
   ZD(2,5) = pars(70);
   ZD(3,0) = pars(71);
   ZD(3,1) = pars(72);
   ZD(3,2) = pars(73);
   ZD(3,3) = pars(74);
   ZD(3,4) = pars(75);
   ZD(3,5) = pars(76);
   ZD(4,0) = pars(77);
   ZD(4,1) = pars(78);
   ZD(4,2) = pars(79);
   ZD(4,3) = pars(80);
   ZD(4,4) = pars(81);
   ZD(4,5) = pars(82);
   ZD(5,0) = pars(83);
   ZD(5,1) = pars(84);
   ZD(5,2) = pars(85);
   ZD(5,3) = pars(86);
   ZD(5,4) = pars(87);
   ZD(5,5) = pars(88);
   ZV(0,0) = pars(89);
   ZV(0,1) = pars(90);
   ZV(0,2) = pars(91);
   ZV(1,0) = pars(92);
   ZV(1,1) = pars(93);
   ZV(1,2) = pars(94);
   ZV(2,0) = pars(95);
   ZV(2,1) = pars(96);
   ZV(2,2) = pars(97);
   ZU(0,0) = pars(98);
   ZU(0,1) = pars(99);
   ZU(0,2) = pars(100);
   ZU(0,3) = pars(101);
   ZU(0,4) = pars(102);
   ZU(0,5) = pars(103);
   ZU(1,0) = pars(104);
   ZU(1,1) = pars(105);
   ZU(1,2) = pars(106);
   ZU(1,3) = pars(107);
   ZU(1,4) = pars(108);
   ZU(1,5) = pars(109);
   ZU(2,0) = pars(110);
   ZU(2,1) = pars(111);
   ZU(2,2) = pars(112);
   ZU(2,3) = pars(113);
   ZU(2,4) = pars(114);
   ZU(2,5) = pars(115);
   ZU(3,0) = pars(116);
   ZU(3,1) = pars(117);
   ZU(3,2) = pars(118);
   ZU(3,3) = pars(119);
   ZU(3,4) = pars(120);
   ZU(3,5) = pars(121);
   ZU(4,0) = pars(122);
   ZU(4,1) = pars(123);
   ZU(4,2) = pars(124);
   ZU(4,3) = pars(125);
   ZU(4,4) = pars(126);
   ZU(4,5) = pars(127);
   ZU(5,0) = pars(128);
   ZU(5,1) = pars(129);
   ZU(5,2) = pars(130);
   ZU(5,3) = pars(131);
   ZU(5,4) = pars(132);
   ZU(5,5) = pars(133);
   ZE(0,0) = pars(134);
   ZE(0,1) = pars(135);
   ZE(0,2) = pars(136);
   ZE(0,3) = pars(137);
   ZE(0,4) = pars(138);
   ZE(0,5) = pars(139);
   ZE(1,0) = pars(140);
   ZE(1,1) = pars(141);
   ZE(1,2) = pars(142);
   ZE(1,3) = pars(143);
   ZE(1,4) = pars(144);
   ZE(1,5) = pars(145);
   ZE(2,0) = pars(146);
   ZE(2,1) = pars(147);
   ZE(2,2) = pars(148);
   ZE(2,3) = pars(149);
   ZE(2,4) = pars(150);
   ZE(2,5) = pars(151);
   ZE(3,0) = pars(152);
   ZE(3,1) = pars(153);
   ZE(3,2) = pars(154);
   ZE(3,3) = pars(155);
   ZE(3,4) = pars(156);
   ZE(3,5) = pars(157);
   ZE(4,0) = pars(158);
   ZE(4,1) = pars(159);
   ZE(4,2) = pars(160);
   ZE(4,3) = pars(161);
   ZE(4,4) = pars(162);
   ZE(4,5) = pars(163);
   ZE(5,0) = pars(164);
   ZE(5,1) = pars(165);
   ZE(5,2) = pars(166);
   ZE(5,3) = pars(167);
   ZE(5,4) = pars(168);
   ZE(5,5) = pars(169);
   ZH(0,0) = pars(170);
   ZH(0,1) = pars(171);
   ZH(0,2) = pars(172);
   ZH(1,0) = pars(173);
   ZH(1,1) = pars(174);
   ZH(1,2) = pars(175);
   ZH(2,0) = pars(176);
   ZH(2,1) = pars(177);
   ZH(2,2) = pars(178);
   ZA(0,0) = pars(179);
   ZA(0,1) = pars(180);
   ZA(0,2) = pars(181);
   ZA(1,0) = pars(182);
   ZA(1,1) = pars(183);
   ZA(1,2) = pars(184);
   ZA(2,0) = pars(185);
   ZA(2,1) = pars(186);
   ZA(2,2) = pars(187);
   ZP(0,0) = pars(188);
   ZP(0,1) = pars(189);
   ZP(1,0) = pars(190);
   ZP(1,1) = pars(191);
   ZN(0,0) = std::complex<double>(pars(192), pars(193));
   ZN(0,1) = std::complex<double>(pars(194), pars(195));
   ZN(0,2) = std::complex<double>(pars(196), pars(197));
   ZN(0,3) = std::complex<double>(pars(198), pars(199));
   ZN(0,4) = std::complex<double>(pars(200), pars(201));
   ZN(1,0) = std::complex<double>(pars(202), pars(203));
   ZN(1,1) = std::complex<double>(pars(204), pars(205));
   ZN(1,2) = std::complex<double>(pars(206), pars(207));
   ZN(1,3) = std::complex<double>(pars(208), pars(209));
   ZN(1,4) = std::complex<double>(pars(210), pars(211));
   ZN(2,0) = std::complex<double>(pars(212), pars(213));
   ZN(2,1) = std::complex<double>(pars(214), pars(215));
   ZN(2,2) = std::complex<double>(pars(216), pars(217));
   ZN(2,3) = std::complex<double>(pars(218), pars(219));
   ZN(2,4) = std::complex<double>(pars(220), pars(221));
   ZN(3,0) = std::complex<double>(pars(222), pars(223));
   ZN(3,1) = std::complex<double>(pars(224), pars(225));
   ZN(3,2) = std::complex<double>(pars(226), pars(227));
   ZN(3,3) = std::complex<double>(pars(228), pars(229));
   ZN(3,4) = std::complex<double>(pars(230), pars(231));
   ZN(4,0) = std::complex<double>(pars(232), pars(233));
   ZN(4,1) = std::complex<double>(pars(234), pars(235));
   ZN(4,2) = std::complex<double>(pars(236), pars(237));
   ZN(4,3) = std::complex<double>(pars(238), pars(239));
   ZN(4,4) = std::complex<double>(pars(240), pars(241));
   UM(0,0) = std::complex<double>(pars(242), pars(243));
   UM(0,1) = std::complex<double>(pars(244), pars(245));
   UM(1,0) = std::complex<double>(pars(246), pars(247));
   UM(1,1) = std::complex<double>(pars(248), pars(249));
   UP(0,0) = std::complex<double>(pars(250), pars(251));
   UP(0,1) = std::complex<double>(pars(252), pars(253));
   UP(1,0) = std::complex<double>(pars(254), pars(255));
   UP(1,1) = std::complex<double>(pars(256), pars(257));
   ZEL(0,0) = std::complex<double>(pars(258), pars(259));
   ZEL(0,1) = std::complex<double>(pars(260), pars(261));
   ZEL(0,2) = std::complex<double>(pars(262), pars(263));
   ZEL(1,0) = std::complex<double>(pars(264), pars(265));
   ZEL(1,1) = std::complex<double>(pars(266), pars(267));
   ZEL(1,2) = std::complex<double>(pars(268), pars(269));
   ZEL(2,0) = std::complex<double>(pars(270), pars(271));
   ZEL(2,1) = std::complex<double>(pars(272), pars(273));
   ZEL(2,2) = std::complex<double>(pars(274), pars(275));
   ZER(0,0) = std::complex<double>(pars(276), pars(277));
   ZER(0,1) = std::complex<double>(pars(278), pars(279));
   ZER(0,2) = std::complex<double>(pars(280), pars(281));
   ZER(1,0) = std::complex<double>(pars(282), pars(283));
   ZER(1,1) = std::complex<double>(pars(284), pars(285));
   ZER(1,2) = std::complex<double>(pars(286), pars(287));
   ZER(2,0) = std::complex<double>(pars(288), pars(289));
   ZER(2,1) = std::complex<double>(pars(290), pars(291));
   ZER(2,2) = std::complex<double>(pars(292), pars(293));
   ZDL(0,0) = std::complex<double>(pars(294), pars(295));
   ZDL(0,1) = std::complex<double>(pars(296), pars(297));
   ZDL(0,2) = std::complex<double>(pars(298), pars(299));
   ZDL(1,0) = std::complex<double>(pars(300), pars(301));
   ZDL(1,1) = std::complex<double>(pars(302), pars(303));
   ZDL(1,2) = std::complex<double>(pars(304), pars(305));
   ZDL(2,0) = std::complex<double>(pars(306), pars(307));
   ZDL(2,1) = std::complex<double>(pars(308), pars(309));
   ZDL(2,2) = std::complex<double>(pars(310), pars(311));
   ZDR(0,0) = std::complex<double>(pars(312), pars(313));
   ZDR(0,1) = std::complex<double>(pars(314), pars(315));
   ZDR(0,2) = std::complex<double>(pars(316), pars(317));
   ZDR(1,0) = std::complex<double>(pars(318), pars(319));
   ZDR(1,1) = std::complex<double>(pars(320), pars(321));
   ZDR(1,2) = std::complex<double>(pars(322), pars(323));
   ZDR(2,0) = std::complex<double>(pars(324), pars(325));
   ZDR(2,1) = std::complex<double>(pars(326), pars(327));
   ZDR(2,2) = std::complex<double>(pars(328), pars(329));
   ZUL(0,0) = std::complex<double>(pars(330), pars(331));
   ZUL(0,1) = std::complex<double>(pars(332), pars(333));
   ZUL(0,2) = std::complex<double>(pars(334), pars(335));
   ZUL(1,0) = std::complex<double>(pars(336), pars(337));
   ZUL(1,1) = std::complex<double>(pars(338), pars(339));
   ZUL(1,2) = std::complex<double>(pars(340), pars(341));
   ZUL(2,0) = std::complex<double>(pars(342), pars(343));
   ZUL(2,1) = std::complex<double>(pars(344), pars(345));
   ZUL(2,2) = std::complex<double>(pars(346), pars(347));
   ZUR(0,0) = std::complex<double>(pars(348), pars(349));
   ZUR(0,1) = std::complex<double>(pars(350), pars(351));
   ZUR(0,2) = std::complex<double>(pars(352), pars(353));
   ZUR(1,0) = std::complex<double>(pars(354), pars(355));
   ZUR(1,1) = std::complex<double>(pars(356), pars(357));
   ZUR(1,2) = std::complex<double>(pars(358), pars(359));
   ZUR(2,0) = std::complex<double>(pars(360), pars(361));
   ZUR(2,1) = std::complex<double>(pars(362), pars(363));
   ZUR(2,2) = std::complex<double>(pars(364), pars(365));
   ZZ(0,0) = pars(366);
   ZZ(0,1) = pars(367);
   ZZ(1,0) = pars(368);
   ZZ(1,1) = pars(369);

}

Eigen::ArrayXd NMSSMtower_physical::get_masses() const
{
   Eigen::ArrayXd pars(53);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MSd(0);
   pars(6) = MSd(1);
   pars(7) = MSd(2);
   pars(8) = MSd(3);
   pars(9) = MSd(4);
   pars(10) = MSd(5);
   pars(11) = MSv(0);
   pars(12) = MSv(1);
   pars(13) = MSv(2);
   pars(14) = MSu(0);
   pars(15) = MSu(1);
   pars(16) = MSu(2);
   pars(17) = MSu(3);
   pars(18) = MSu(4);
   pars(19) = MSu(5);
   pars(20) = MSe(0);
   pars(21) = MSe(1);
   pars(22) = MSe(2);
   pars(23) = MSe(3);
   pars(24) = MSe(4);
   pars(25) = MSe(5);
   pars(26) = Mhh(0);
   pars(27) = Mhh(1);
   pars(28) = Mhh(2);
   pars(29) = MAh(0);
   pars(30) = MAh(1);
   pars(31) = MAh(2);
   pars(32) = MHpm(0);
   pars(33) = MHpm(1);
   pars(34) = MChi(0);
   pars(35) = MChi(1);
   pars(36) = MChi(2);
   pars(37) = MChi(3);
   pars(38) = MChi(4);
   pars(39) = MCha(0);
   pars(40) = MCha(1);
   pars(41) = MFe(0);
   pars(42) = MFe(1);
   pars(43) = MFe(2);
   pars(44) = MFd(0);
   pars(45) = MFd(1);
   pars(46) = MFd(2);
   pars(47) = MFu(0);
   pars(48) = MFu(1);
   pars(49) = MFu(2);
   pars(50) = MVWm;
   pars(51) = MVP;
   pars(52) = MVZ;

   return pars;
}

void NMSSMtower_physical::set_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MFv(0) = pars(2);
   MFv(1) = pars(3);
   MFv(2) = pars(4);
   MSd(0) = pars(5);
   MSd(1) = pars(6);
   MSd(2) = pars(7);
   MSd(3) = pars(8);
   MSd(4) = pars(9);
   MSd(5) = pars(10);
   MSv(0) = pars(11);
   MSv(1) = pars(12);
   MSv(2) = pars(13);
   MSu(0) = pars(14);
   MSu(1) = pars(15);
   MSu(2) = pars(16);
   MSu(3) = pars(17);
   MSu(4) = pars(18);
   MSu(5) = pars(19);
   MSe(0) = pars(20);
   MSe(1) = pars(21);
   MSe(2) = pars(22);
   MSe(3) = pars(23);
   MSe(4) = pars(24);
   MSe(5) = pars(25);
   Mhh(0) = pars(26);
   Mhh(1) = pars(27);
   Mhh(2) = pars(28);
   MAh(0) = pars(29);
   MAh(1) = pars(30);
   MAh(2) = pars(31);
   MHpm(0) = pars(32);
   MHpm(1) = pars(33);
   MChi(0) = pars(34);
   MChi(1) = pars(35);
   MChi(2) = pars(36);
   MChi(3) = pars(37);
   MChi(4) = pars(38);
   MCha(0) = pars(39);
   MCha(1) = pars(40);
   MFe(0) = pars(41);
   MFe(1) = pars(42);
   MFe(2) = pars(43);
   MFd(0) = pars(44);
   MFd(1) = pars(45);
   MFd(2) = pars(46);
   MFu(0) = pars(47);
   MFu(1) = pars(48);
   MFu(2) = pars(49);
   MVWm = pars(50);
   MVP = pars(51);
   MVZ = pars(52);

}

void NMSSMtower_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';

   ostr << "----------------------------------------\n"
           "pole mass mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZD = " << ZD << '\n';
   ostr << "ZV = " << ZV << '\n';
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';
   ostr << "ZEL = " << ZEL << '\n';
   ostr << "ZER = " << ZER << '\n';
   ostr << "ZDL = " << ZDL << '\n';
   ostr << "ZDR = " << ZDR << '\n';
   ostr << "ZUL = " << ZUL << '\n';
   ostr << "ZUR = " << ZUR << '\n';
   ostr << "ZZ = " << ZZ << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const NMSSMtower_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
