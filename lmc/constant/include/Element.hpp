/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 9/23/20 1:29 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/27/23 9:50 PM                                                           *
 **************************************************************************************************/

/*! \file  Element.hpp
 *  \brief File for the Element class definition and implementation.
 */

#ifndef LMC_CONSTANT_INCLUDE_ELEMENT_HPP_
#define LMC_CONSTANT_INCLUDE_ELEMENT_HPP_

#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
#include <boost/functional/hash.hpp>


/*! \enum ElementName
 *  \brief Define the names of the elements that can be used. X is used for the vacancy.
 */
enum class ElementName
{
  X,
  H,
  He,
  Li,
  Be,
  B,
  C,
  N,
  O,
  F,
  Ne,
  Na,
  Mg,
  Al,
  Si,
  P,
  S,
  Cl,
  Ar,
  K,
  Ca,
  Sc,
  Ti,
  V,
  Cr,
  Mn,
  Fe,
  Co,
  Ni,
  Cu,
  Zn,
  Ga,
  Ge,
  As,
  Se,
  Br,
  Kr,
  Rb,
  Sr,
  Y,
  Zr,
  Nb,
  Mo,
  Tc,
  Ru,
  Rh,
  Pd,
  Ag,
  Cd,
  In,
  Sn,
  Sb,
  Te,
  I,
  Xe,
  Cs,
  Ba,
  La,
  Ce,
  Pr,
  Nd,
  Pm,
  Sm,
  Eu,
  Gd,
  Tb,
  Dy,
  Ho,
  Er,
  Tm,
  Yb,
  Lu,
  Hf,
  Ta,
  W,
  Re,
  Os,
  Ir,
  Pt,
  Au,
  Hg,
  Tl,
  Pb,
  Bi,
  Po,
  At,
  Rn,
  Fr,
  Ra,
  Ac,
  Th,
  Pa,
  U,
  Np,
  Pu,
  Am,
  Cm,
  Bk,
  Cf,
  Es,
  Fm,
  Md,
  No,
  Lr,
  Rf,
  Db,
  Sg,
  Bh,
  Hs,
  Mt,
  Ds,
  Rg,
  Cn,
  Nh,
  Fl,
  Mc,
  Lv,
  Ts,
  Og,
};

/*! \class Element
 *  \brief A data structure to represent a chemical element.
 *
 *  This class encapsulates the data and operations for a chemical element.
 */
class Element
{
public:
  /*! \brief Default constructor.
   */
  Element() = default;

  /*! \brief Constructor to create an Element instance from an ElementName.
   *  \param element_type : The symbol of this element.
   */
  constexpr explicit Element(ElementName element_type) : element_name_(element_type) {}

  /*! \brief Constructor to create an Element instance from a string.
   *  \param element_string : The name of this element.
   */
  explicit Element(const std::string &element_string)
  {
    auto it = element_string_map_.find(element_string);
    // If the string does not correspond to a valid element, throw an exception
    if (it == element_string_map_.end())
    {
      throw std::invalid_argument("Invalid element string: " + element_string);
    }
    element_name_ = it->second;
  }

  /*! \brief Comparison and other operator overloads.
   */
  constexpr explicit operator ElementName() const { return element_name_; }
  explicit operator bool() = delete;
  constexpr bool operator==(Element rhs) const
  {
    return element_name_ == rhs.element_name_;
  }
  constexpr bool operator!=(Element rhs) const
  {
    return element_name_ != rhs.element_name_;
  }
  constexpr bool operator==(ElementName rhs) const
  {
    return element_name_ == rhs;
  }
  constexpr bool operator!=(ElementName rhs) const
  {
    return element_name_ != rhs;
  }
  bool operator<(Element rhs) const
  {
    return GetElementString() < rhs.GetElementString();
  }
  bool operator<(const ElementName rhs) const
  {
    return GetElementString() < Element(rhs).GetElementString();
  }

  
  /*! \brief Hash function.
   */
  friend size_t hash_value(Element element)
  {
    return static_cast<std::size_t>(element.element_name_);
  }

  /*! \brief Returns the atomic index of this element.
   *  \return : The atomic index of this element.
   */
  [[nodiscard]] size_t GetAtomicIndex() const
  {
    return GetProperties().index;
  }

  /*! \brief Returns the name of this element.
   *  \return : The name of this element.
   */
  [[nodiscard]] std::string GetElementString() const
  {
    return GetProperties().name;
  }

  /*! \brief Returns the mass of this element.
   *  \return : The mass of this element.
   */
  [[nodiscard]] double GetMass() const
  {
    return GetProperties().mass;
  }

  /*! \brief Returns the electronegativity of this element.
   *  \return : The electronegativity of this element.
   */
  [[nodiscard]] double GetElectronegativity() const
  {
    return GetProperties().electronegativity;
  }

  /*! \brief Returns the valence number of electron of element of this element.
   *  \return : The valence number of electron (Sv) of this element.
   */
  [[nodiscard]] size_t GetSv() const
  {

    auto it = element_Sv_map_.find(element_name_);
    if (it == element_Sv_map_.end())
    {
      throw std::invalid_argument("Element not found in the Sv Map");
    }
    return it->second;
  }

  /*! \brief stream operator output for this element.
   *  \param os      : The output stream.
   *  \param element : The Element to be streamed.
   *  \return        : The output stream.
   */
  friend std::ostream &operator<<(std::ostream &os, const Element &element)
  {
    os << element.GetElementString();
    return os;
  }

private:
  /*! \struct ElementProperties
   *  \brief A data structure to hold the properties of an element.
   */
  struct ElementProperties
  {
    /// The atomic index of this element
    size_t index;
    /// The symbol of this element in string.
    std::string name;
    /// The mass of this element.
    double mass;
    /// Electronegativity of this element on Pauling Scale.
    /// 0 for no value.
    /// https://pubchem.ncbi.nlm.nih.gov/ptable/electronegativity/
    double electronegativity;
  };

  /// Inline static map holding the properties of each ElementName.
  inline static const std::unordered_map<ElementName, ElementProperties> element_properties_map_ = {
      {ElementName::X, {0, "X", 0.0, 1.0}},
      {ElementName::H, {1, "H", 1.008, 2.2}},
      {ElementName::He, {2, "He", 4.0026, 0.0}},
      {ElementName::Li, {3, "Li", 6.94, 0.98}},
      {ElementName::Be, {4, "Be", 9.0122, 1.57}},
      {ElementName::B, {5, "B", 10.81, 2.04}},
      {ElementName::C, {6, "C", 12.011, 2.55}},
      {ElementName::N, {7, "N", 14.007, 3.04}},
      {ElementName::O, {8, "O", 15.999, 3.44}},
      {ElementName::F, {9, "F", 18.998, 3.98}},
      {ElementName::Ne, {10, "Ne", 20.180, 0.0}},
      {ElementName::Na, {11, "Na", 22.990, 0.93}},
      {ElementName::Mg, {12, "Mg", 24.305, 1.31}},
      {ElementName::Al, {13, "Al", 26.982, 1.61}},
      {ElementName::Si, {14, "Si", 28.085, 1.9}},
      {ElementName::P, {15, "P", 30.974, 2.19}},
      {ElementName::S, {16, "S", 32.06, 2.58}},
      {ElementName::Cl, {17, "Cl", 35.45, 3.16}},
      {ElementName::Ar, {18, "Ar", 39.948, 0.0}},
      {ElementName::K, {19, "K", 39.098, 0.82}},
      {ElementName::Ca, {20, "Ca", 40.078, 1.0}},
      {ElementName::Sc, {21, "Sc", 44.956, 1.36}},
      {ElementName::Ti, {22, "Ti", 47.867, 1.54}},
      {ElementName::V, {23, "V", 50.942, 1.63}},
      {ElementName::Cr, {24, "Cr", 51.996, 1.66}},
      {ElementName::Mn, {25, "Mn", 54.938, 1.55}},
      {ElementName::Fe, {26, "Fe", 55.845, 1.83}},
      {ElementName::Co, {27, "Co", 58.933, 1.88}},
      {ElementName::Ni, {28, "Ni", 58.693, 1.91}},
      {ElementName::Cu, {29, "Cu", 63.546, 1.9}},
      {ElementName::Zn, {30, "Zn", 65.38, 1.65}},
      {ElementName::Ga, {31, "Ga", 69.723, 1.81}},
      {ElementName::Ge, {32, "Ge", 72.630, 2.01}},
      {ElementName::As, {33, "As", 74.922, 2.18}},
      {ElementName::Se, {34, "Se", 78.971, 2.55}},
      {ElementName::Br, {35, "Br", 79.904, 2.96}},
      {ElementName::Kr, {36, "Kr", 83.798, 3.0}},
      {ElementName::Rb, {37, "Rb", 85.468, 0.82}},
      {ElementName::Sr, {38, "Sr", 87.62, 0.95}},
      {ElementName::Y, {39, "Y", 88.906, 1.22}},
      {ElementName::Zr, {40, "Zr", 91.224, 1.33}},
      {ElementName::Nb, {41, "Nb", 92.906, 1.6}},
      {ElementName::Mo, {42, "Mo", 95.95, 2.16}},
      {ElementName::Tc, {43, "Tc", 98, 1.9}},
      {ElementName::Ru, {44, "Ru", 101.07, 2.2}},
      {ElementName::Rh, {45, "Rh", 102.91, 2.28}},
      {ElementName::Pd, {46, "Pd", 106.42, 2.2}},
      {ElementName::Ag, {47, "Ag", 107.87, 1.93}},
      {ElementName::Cd, {48, "Cd", 112.41, 1.69}},
      {ElementName::In, {49, "In", 114.82, 1.78}},
      {ElementName::Sn, {50, "Sn", 118.71, 1.96}},
      {ElementName::Sb, {51, "Sb", 121.76, 2.05}},
      {ElementName::Te, {52, "Te", 127.60, 2.1}},
      {ElementName::I, {53, "I", 126.90, 2.66}},
      {ElementName::Xe, {54, "Xe", 131.29, 2.6}},
      {ElementName::Cs, {55, "Cs", 132.91, 0.79}},
      {ElementName::Ba, {56, "Ba", 137.33, 0.89}},
      {ElementName::La, {57, "La", 138.91, 1.1}},
      {ElementName::Ce, {58, "Ce", 140.12, 1.12}},
      {ElementName::Pr, {59, "Pr", 140.91, 1.13}},
      {ElementName::Nd, {60, "Nd", 144.24, 1.14}},
      {ElementName::Pm, {61, "Pm", 145, 0.0}},
      {ElementName::Sm, {62, "Sm", 150.36, 1.17}},
      {ElementName::Eu, {63, "Eu", 151.96, 0.0}},
      {ElementName::Gd, {64, "Gd", 157.25, 1.2}},
      {ElementName::Tb, {65, "Tb", 158.93, 0.0}},
      {ElementName::Dy, {66, "Dy", 162.50, 1.22}},
      {ElementName::Ho, {67, "Ho", 164.93, 1.23}},
      {ElementName::Er, {68, "Er", 167.26, 1.24}},
      {ElementName::Tm, {69, "Tm", 168.93, 1.25}},
      {ElementName::Yb, {70, "Yb", 173.05, 0.0}},
      {ElementName::Lu, {71, "Lu", 174.97, 1.27}},
      {ElementName::Hf, {72, "Hf", 178.49, 1.3}},
      {ElementName::Ta, {73, "Ta", 180.95, 1.5}},
      {ElementName::W, {74, "W", 183.84, 2.36}},
      {ElementName::Re, {75, "Re", 186.21, 1.9}},
      {ElementName::Os, {76, "Os", 190.23, 2.2}},
      {ElementName::Ir, {77, "Ir", 192.22, 2.2}},
      {ElementName::Pt, {78, "Pt", 195.08, 2.28}},
      {ElementName::Au, {79, "Au", 196.97, 2.54}},
      {ElementName::Hg, {80, "Hg", 200.59, 2}},
      {ElementName::Tl, {81, "Tl", 204.38, 1.62}},
      {ElementName::Pb, {82, "Pb", 207.2, 2.33}},
      {ElementName::Bi, {83, "Bi", 208.98, 2.02}},
      {ElementName::Po, {84, "Po", 209, 2.0}},
      {ElementName::At, {85, "At", 210, 2.2}},
      {ElementName::Rn, {86, "Rn", 222, 0.0}},
      {ElementName::Fr, {87, "Fr", 223, 0.7}},
      {ElementName::Ra, {88, "Ra", 226, 0.9}},
      {ElementName::Ac, {89, "Ac", 227, 1.1}},
      {ElementName::Th, {90, "Th", 232.04, 1.3}},
      {ElementName::Pa, {91, "Pa", 231.04, 1.5}},
      {ElementName::U, {92, "U", 238.03, 1.38}},
      {ElementName::Np, {93, "Np", 237, 1.36}},
      {ElementName::Pu, {94, "Pu", 244, 1.28}},
      {ElementName::Am, {95, "Am", 243, 1.3}},
      {ElementName::Cm, {96, "Cm", 247, 1.3}},
      {ElementName::Bk, {97, "Bk", 247, 1.3}},
      {ElementName::Cf, {98, "Cf", 251, 1.3}},
      {ElementName::Es, {99, "Es", 252, 1.3}},
      {ElementName::Fm, {100, "Fm", 257, 1.3}},
      {ElementName::Md, {101, "Md", 258, 1.3}},
      {ElementName::No, {102, "No", 259, 1.3}},
      {ElementName::Lr, {103, "Lr", 266, 1.3}},
      {ElementName::Rf, {104, "Rf", 267, 0.0}},
      {ElementName::Db, {105, "Db", 270, 0.0}},
      {ElementName::Sg, {106, "Sg", 271, 0.0}},
      {ElementName::Bh, {107, "Bh", 270, 0.0}},
      {ElementName::Hs, {108, "Hs", 277, 0.0}},
      {ElementName::Mt, {109, "Mt", 276, 0.0}},
      {ElementName::Ds, {110, "Ds", 281, 0.0}},
      {ElementName::Rg, {111, "Rg", 282, 0.0}},
      {ElementName::Cn, {112, "Cn", 285, 0.0}},
      {ElementName::Nh, {113, "Nh", 286, 0.0}},
      {ElementName::Fl, {114, "Fl", 289, 0.0}},
      {ElementName::Mc, {115, "Mc", 288, 0.0}},
      {ElementName::Lv, {116, "Lv", 293, 0.0}},
      {ElementName::Ts, {117, "Ts", 294, 0.0}},
      {ElementName::Og, {118, "Og", 294, 0.0}},
  };

  /// Inline static map to map string names to ElementName.
  inline static const std::unordered_map<std::string, ElementName> element_string_map_ = {
      {"X", ElementName::X},
      {"H", ElementName::H},
      {"He", ElementName::He},
      {"Li", ElementName::Li},
      {"Be", ElementName::Be},
      {"B", ElementName::B},
      {"C", ElementName::C},
      {"N", ElementName::N},
      {"O", ElementName::O},
      {"F", ElementName::F},
      {"Ne", ElementName::Ne},
      {"Na", ElementName::Na},
      {"Mg", ElementName::Mg},
      {"Al", ElementName::Al},
      {"Si", ElementName::Si},
      {"P", ElementName::P},
      {"S", ElementName::S},
      {"Cl", ElementName::Cl},
      {"Ar", ElementName::Ar},
      {"K", ElementName::K},
      {"Ca", ElementName::Ca},
      {"Sc", ElementName::Sc},
      {"Ti", ElementName::Ti},
      {"V", ElementName::V},
      {"Cr", ElementName::Cr},
      {"Mn", ElementName::Mn},
      {"Fe", ElementName::Fe},
      {"Co", ElementName::Co},
      {"Ni", ElementName::Ni},
      {"Cu", ElementName::Cu},
      {"Zn", ElementName::Zn},
      {"Ga", ElementName::Ga},
      {"Ge", ElementName::Ge},
      {"As", ElementName::As},
      {"Se", ElementName::Se},
      {"Br", ElementName::Br},
      {"Kr", ElementName::Kr},
      {"Rb", ElementName::Rb},
      {"Sr", ElementName::Sr},
      {"Y", ElementName::Y},
      {"Zr", ElementName::Zr},
      {"Nb", ElementName::Nb},
      {"Mo", ElementName::Mo},
      {"Tc", ElementName::Tc},
      {"Ru", ElementName::Ru},
      {"Rh", ElementName::Rh},
      {"Pd", ElementName::Pd},
      {"Ag", ElementName::Ag},
      {"Cd", ElementName::Cd},
      {"In", ElementName::In},
      {"Sn", ElementName::Sn},
      {"Sb", ElementName::Sb},
      {"Te", ElementName::Te},
      {"I", ElementName::I},
      {"Xe", ElementName::Xe},
      {"Cs", ElementName::Cs},
      {"Ba", ElementName::Ba},
      {"La", ElementName::La},
      {"Ce", ElementName::Ce},
      {"Pr", ElementName::Pr},
      {"Nd", ElementName::Nd},
      {"Pm", ElementName::Pm},
      {"Sm", ElementName::Sm},
      {"Eu", ElementName::Eu},
      {"Gd", ElementName::Gd},
      {"Tb", ElementName::Tb},
      {"Dy", ElementName::Dy},
      {"Ho", ElementName::Ho},
      {"Er", ElementName::Er},
      {"Tm", ElementName::Tm},
      {"Yb", ElementName::Yb},
      {"Lu", ElementName::Lu},
      {"Hf", ElementName::Hf},
      {"Ta", ElementName::Ta},
      {"W", ElementName::W},
      {"Re", ElementName::Re},
      {"Os", ElementName::Os},
      {"Ir", ElementName::Ir},
      {"Pt", ElementName::Pt},
      {"Au", ElementName::Au},
      {"Hg", ElementName::Hg},
      {"Tl", ElementName::Tl},
      {"Pb", ElementName::Pb},
      {"Bi", ElementName::Bi},
      {"Po", ElementName::Po},
      {"At", ElementName::At},
      {"Rn", ElementName::Rn},
      {"Fr", ElementName::Fr},
      {"Ra", ElementName::Ra},
      {"Ac", ElementName::Ac},
      {"Th", ElementName::Th},
      {"Pa", ElementName::Pa},
      {"U", ElementName::U},
      {"Np", ElementName::Np},
      {"Pu", ElementName::Pu},
      {"Am", ElementName::Am},
      {"Cm", ElementName::Cm},
      {"Bk", ElementName::Bk},
      {"Cf", ElementName::Cf},
      {"Es", ElementName::Es},
      {"Fm", ElementName::Fm},
      {"Md", ElementName::Md},
      {"No", ElementName::No},
      {"Lr", ElementName::Lr},
      {"Rf", ElementName::Rf},
      {"Db", ElementName::Db},
      {"Sg", ElementName::Sg},
      {"Bh", ElementName::Bh},
      {"Hs", ElementName::Hs},
      {"Mt", ElementName::Mt},
      {"Ds", ElementName::Ds},
      {"Rg", ElementName::Rg},
      {"Cn", ElementName::Cn},
      {"Nh", ElementName::Nh},
      {"Fl", ElementName::Fl},
      {"Mc", ElementName::Mc},
      {"Lv", ElementName::Lv},
      {"Ts", ElementName::Ts},
      {"Og", ElementName::Og},
  };

  /// Inline static map holding the valence number of few elements.
  inline static const std::unordered_map<ElementName, size_t> element_Sv_map_ = {
      {ElementName::X, 1},
      {ElementName::Al, 3},
      {ElementName::Ti, 4},
      {ElementName::V, 5},
      {ElementName::Ni, 10},
      {ElementName::Zr, 4},
      {ElementName::Nb, 5},
      {ElementName::Mo, 6},
      {ElementName::Hf, 4},
      {ElementName::Ta, 5},
      {ElementName::W, 6},

  };

  /*! \brief Returns the properties of this element.
   *  \return : The properties of this element.
   */
  [[nodiscard]] const ElementProperties &GetProperties() const
  {
    auto it = element_properties_map_.find(element_name_);
    if (it == element_properties_map_.end())
    {
      throw std::invalid_argument("Unexpected element");
    }
    return it->second;
  }

  /// The symbol of this element.
  ElementName element_name_{};
};

#endif // LMC_CONSTANT_INCLUDE_ELEMENT_HPP_
