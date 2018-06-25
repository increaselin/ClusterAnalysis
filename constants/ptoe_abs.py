#! /usr/bin/env python
#
# 
# Module: ptoe_abs.py
#
# Creator: Mike Lum
#
# Description: Abundance and look-up table for elements.
#
# Contents: 
#   (Data and LUT only)
#
# Revision History:
#    Date        Vers.    Author        Description
#    4/05/2018   1.0f1    Lum         Formatting and cleanup for github
#
# To Do:
#    - This should be in a real (sqlite3) DB...
#


# If you _REALLY_ need an ionization state greater than 15, feel
# free to enter it here:
RomanNumberLUT = {'I':1,'II':2,'III':3,'IV':4,'V':5,'VI':6,'VII':7,'VIII':8,\
                 'IX':9,'X':10,'XI':11,'XII':12,'XIII':13,'XIV':14,'XV':15}
RomanNumbers = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII',\
                'XIII','XIV','XV']

# Periodic Table format: 
#'At. No.','Symbol','Name','Ab. Solar','Ab. Arcturus'
# Ab. Solar from: (Asplund',2009)
# Ab. Arcturus from: (Ramirez & Allende-Prieto, 2011)
# Abundance adjustments for Arcturus are only available for:
# C(6), O(8), Na(11), Mg(12), Al(13), Si(14), K(19), Ca(20), Sc(21)*, Ti(22)*,
# V(23), Cr(24), Mn(25), Fe(26), Co(27), Ni(28) and Zn(30). 
# We have used our own measurements for Arcturus for Ba(56).
# All other Arcturus "adjustments" are taken from the model (Table 2) used by
# Peterson, Dalle Ore, and Kurucz (1993). 
# *-For Sc, we adopt the Ramirez et al. value for ScII, since the ionized 
# species has been a better indicator in our previous work. 
# Similarly for the TiI value.
# Note that you can sorta assume that the element entries are in order, and
# that the index (-1, since Python arrays are 0-based) is also the At. No.

# Indices
atnoIdx = 0
symbolIdx = 1
nameIdx = 2
abIdx = 3
solarAbIdx = 3
solarErrIdx = 4
giantAbIdx = 5

ptoe = [
(1,  'H', 'Hydrogen',      12.00, 0.00,  12.00),
(2,  'He','Helium',        10.93, 0.01,  10.93),
(3,  'Li','Lithium',        3.26, 0.10,   0.66),
(4,  'Be','Beryllium',      1.38, 0.09,   0.65),
(5,  'B', 'Boron',          2.70, 0.20,   2.10),
(6,  'C', 'Carbon',         8.43, 0.05,   8.34),
(7,  'N', 'Nitrogen',       7.83, 0.05,   7.85),
(8,  'O', 'Oxygen',         8.69, 0.05,   8.67),
(9,  'F', 'Fluorine',       4.56, 0.30,   4.06),
(10, 'Ne','Neon',           7.93, 0.10,   7.59),
(11, 'Na','Sodium',         6.24, 0.04,   5.83),
(12, 'Mg','Magnesium',      7.60, 0.04,   7.45),
(13, 'Al','Aluminium',      6.45, 0.03,   6.28),
(14, 'Si','Silicon',        7.51, 0.03,   7.32),
(15, 'P', 'Phosphorus',     5.41, 0.03,   4.95),
(16, 'S', 'Sulphur',        7.12, 0.03,   7.26),
(17, 'Cl','Chlorine',       5.50, 0.30,   5.00),
(18, 'Ar','Argon',          6.40, 0.13,   6.06),
(19, 'K', 'Potassium',      5.03, 0.03,   4.71),
(20, 'Ca','Calcium',        6.34, 0.04,   5.93),
(21, 'Sc','Scandium',       3.15, 0.04,   2.86),
(22, 'Ti','Titanium',       4.95, 0.05,   4.70),
(23, 'V', 'Vanadium',       3.93, 0.08,   3.61),
(24, 'Cr','Chromium',       5.64, 0.04,   5.07),
(25, 'Mn','Manganese',      5.43, 0.05,   4.70),
(26, 'Fe','Iron',           7.50, 0.04,   6.98),
(27, 'Co','Cobalt',         4.99, 0.07,   4.56),
(28, 'Ni','Nickel',         6.22, 0.04,   5.76),
(29, 'Cu','Copper',         4.19, 0.04,   3.71),
(30, 'Zn','Zinc',           4.56, 0.05,   4.26),
(31, 'Ga','Gallium',        3.04, 0.09,   2.38),
(32, 'Ge','Germanium',      3.65, 0.10,   2.91),
(33, 'As','Arsenic',       -9.99,-9.99,   1.87),
(34, 'Se','Selenium',      -9.99,-9.99,   2.85),
(35, 'Br','Bromine',       -9.99,-9.99,   2.13),
(36, 'Kr','Krypton',        3.25, 0.06,   2.73),
(37, 'Rb','Rubidium',       2.52, 0.10,   2.10),
(38, 'Sr','Strontium',      2.87, 0.07,   2.40),
(39, 'Y', 'Yttrium',        2.21, 0.05,   1.74),
(40, 'Zr','Zirconium',      2.58, 0.04,   2.00),
(41, 'Nb','Niobium',        1.46, 0.04,   0.92),
(42, 'Mo','Molybdenum',     1.88, 0.08,   1.42),
(43, 'Tc','Technetium',    -9.99,-9.99,  -9.99),
(44, 'Ru','Ruthenium',      1.75, 0.08,   1.34),
(45, 'Rh','Rhodium',        0.91, 0.10,   0.62),
(46, 'Pd','Palladium',      1.57, 0.10,   1.19),
(47, 'Ag','Silver',         0.94, 0.10,   0.44),
(48, 'Cd','Cadmium',       -9.99,-9.99,   1.36),
(49, 'In','Indium',         0.80, 0.20,   0.96),
(50, 'Sn','Tin',            2.04, 0.10,   1.50),
(51, 'Sb','Antimony',      -9.99,-9.99,   0.50),
(52, 'Te','Tellurium',     -9.99,-9.99,   1.74),
(53, 'I', 'Iodine',        -9.99,-9.99,   1.01),
(54, 'Xe','Xenon',          2.24, 0.06,   1.73),
(55, 'Cs','Caesium',       -9.99,-9.99,   0.62),
(56, 'Ba','Barium',         2.18, 0.09,   1.63),
(57, 'La','Lanthanum',      1.10, 0.04,   0.72),
(58, 'Ce','Cerium',         1.58, 0.04,   1.05),
(59, 'Pr','Praseodymium',   0.72, 0.04,   0.21),
(60, 'Nd','Neodymium',      1.42, 0.04,   1.00),
(61, 'Pm','Promethium',    -9.99,-9.99,  -9.99),
(62, 'Sm','Samarium',       0.96, 0.04,   0.50),
(63, 'Eu','Europium',       0.52, 0.04,   0.01),
(64, 'Gd','Gadolinium',     1.07, 0.04,   0.62),
(65, 'Tb','Terbium',        0.30, 0.10,  -0.40),
(66, 'Dy','Dysprosium',     1.10, 0.04,   0.60),
(67, 'Ho','Holmium',        0.48, 0.11,  -0.24),
(68, 'Er','Erbium',         0.92, 0.05,   0.43),
(69, 'Tm','Thulium',        0.10, 0.04,  -0.50),
(70, 'Yb','Ytterbium',      0.84, 0.11,   0.58),
(71, 'Lu','Lutetium',       0.10, 0.09,   0.26),
(72, 'Hf','Hafnium',        0.85, 0.04,   0.38),
(73, 'Ta','Tantalium',     -9.99,-9.99,  -0.37),
(74, 'W', 'Tungsten',       0.85, 0.12,   0.61),
(75, 'Re','Rhenium',       -9.99,-9.99,  -0.23),
(76, 'Os','Osmium',         1.40, 0.08,   0.95),
(77, 'Ir','Iridium',        1.38, 0.07,   0.85),
(78, 'Pt','Platinum',      -9.99,-9.99,   1.30),
(79, 'Au','Gold',           0.92, 0.10,   0.51),
(80, 'Hg','Mercury',       -9.99,-9.99,   0.59),
(81, 'Tl','Thallium',       0.90, 0.20,   0.40),
(82, 'Pb','Lead',           1.75, 0.10,   1.35),
(83, 'Bi','Bismuth',       -9.99,-9.99,   0.21),
(84, 'Po','Polonium',      -9.99,-9.99,  -9.99),
(85, 'A','Astatine',       -9.99,-9.99,  -9.99),
(86, 'Rn','Radon',         -9.99,-9.99,  -9.99),
(87, 'Fr','Francium',      -9.99,-9.99,  -9.99),
(88, 'Ra','Radium',        -9.99,-9.99,  -9.99),
(89, 'Ac','Actinium',      -9.99,-9.99,  -9.99),
(90, 'Th','Thorium',        0.02, 0.10,  -0.38),
(91, 'Pa','Protactinium',  -9.99,-9.99,  -9.99),
(92, 'U', 'Uranium',       -9.99,-9.99,  -0.97),
(93, 'Np','Neptunium',     -9.99,-9.99,  -9.99),
(94, 'Pu','Plutonium',     -9.99,-9.99,  -9.99),
(95, 'Am','Americium',     -9.99,-9.99,  -9.99),
(96, 'Cm','Curium',        -9.99,-9.99,  -9.99),
(97, 'Bk','Berkelium',     -9.99,-9.99,  -9.99),
(98, 'Cf','Californium',   -9.99,-9.99,  -9.99),
(99, 'Es','Einsteinium',   -9.99,-9.99,  -9.99),
(100,'Fm','Fermium',       -9.99,-9.99,  -9.99),
(101,'Md','Mendelevium',   -9.99,-9.99,  -9.99),
(102,'No','Nobelium',      -9.99,-9.99,  -9.99),
(103,'Lr','Lawrencium',    -9.99,-9.99,  -9.99),
(104,'Rf','Rutherfordium', -9.99,-9.99,  -9.99),
(105,'Db','Dubnium',       -9.99,-9.99,  -9.99),
(106,'Sg','Seaborgium',    -9.99,-9.99,  -9.99),
(107,'Bh','Bohrium',       -9.99,-9.99,  -9.99),
(108,'Hs','Hassium',       -9.99,-9.99,  -9.99),
(109,'Mt','Meitnerium',    -9.99,-9.99,  -9.99),
(110,'Ds','Darmstadtium',  -9.99,-9.99,  -9.99),
(111,'Rg','Roentgenium',   -9.99,-9.99,  -9.99),
(112,'Cn','Coppernicium',  -9.99,-9.99,  -9.99),
(113,'Nh','Nihonium',      -9.99,-9.99,  -9.99),
(114,'Fl','Flerovium',     -9.99,-9.99,  -9.99),
(115,'Mc','Moscovium',     -9.99,-9.99,  -9.99),
(116,'Lv','Livermorium',   -9.99,-9.99,  -9.99),
(117,'Ts','Tennessine',    -9.99,-9.99,  -9.99),
(118,'Og','Oganesson',     -9.99,-9.99,  -9.99),
(74, 'Wl','Tungsten',       0.85,0.85),]    # Note: this symbol is wrong, due to MOOG being stupid
