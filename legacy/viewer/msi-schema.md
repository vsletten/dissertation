## Cerius2 Molecular Modeling Software

Cerius2 a molecular modeling and simulation software package created by Molecular Simulations Inc. (MSI), a company focused on computational chemistry tools. In the late 1990s, Cerius2 was widely used in academic and industrial settings, including by researchers at institutions like Carnegie Mellon, for studying the properties of chemicals and materials at the molecular level.

The software provided a comprehensive suite of tools for tasks such as molecular visualization, structure building, and forcefield-based simulations. It featured a graphical user interface (GUI) called the Cerius2 Visualizer, which served as the core platform for managing and presenting molecular models. This interface allowed users to build, visualize, and manipulate atomistic models in 3D, supported by various specialized modules for applications like crystallization, polymer modeling, quantum mechanics, and classical simulations. These modules plugged into the Visualizer, offering a consistent, mouse-driven control system.

###  MSI CERIUS2 DataModel File Version 2 0
Below is an example of a Cerius2 DataModel file, which is a text-based format used to store molecular models and simulation data. The file is structured as a tree of objects, with each object containing a set of attributes.
```
# MSI CERIUS2 DataModel File Version 2 0
(1 Model [1 is the object ID, Model is the object type and delimits model data]
 (A C Label methane) [A = attribute tag, C = string-type attribute]
 (2 Atom [2 is the object ID, Atom is the object type and delimits data for first atom]
  (A I ACL "1 H") [A = attribute tag, I = integer-type attribute]
  (A F Charge 0.028) [F = floating-point type attrubute, Charge is the nameof the attribute, 
		the following number is the attribute's value]
  (A C Label H1)
 )
 (3 Atom
  (A I ACL "6 C")
  (A F Charge -0.11)
  (A D XYZ (1.087 0 0)) [D = double-precision type attribute, the three values are an array]
  (A C Label C2)
 )
 (4 Atom
  (A I ACL "1 H")
  (A F Charge 0.028)
  (A D XYZ (1.4493 1.02483 0))
  (A C Label H3)
 )
 (5 Atom
  (A I ACL "1 H")
  (A F Charge 0.028)
  (A D XYZ (1.44936 -0.51245 -0.88749))
  (A C Label H4)
 )
 (6 Atom
  (A I ACL "1 H")
  (A F Charge 0.028)
  (A D XYZ (1.44934 -0.51236 0.88756))
  (A C Label H5)
 )
 (7 Bond [Bond is the object type and delimits data for first bond]
  (A O Atom1 2) [O = object-OD type attribute]
  (A O Atom2 3)
 )
 (8 Bond
  (A O Atom1 3)
  (A O Atom2 4)
 )
 (9 Bond
  (A O Atom1 3)
  (A O Atom2 5)
 )
 (10 Bond
  (A O Atom1 3)
  (A O Atom2 6)
 )
)
```
### Attribute Types and Object Types
Allowable attribute types are: B (byte), C (string), D (double-precision number), F (floating-point number), I (integer), O (object ID), S (short), T (table).

### Object Types
Allowable object types include ACL, Atom, Bond, Group, SCL, Sequence, Subunit, string.

### Attribute Values
Allowable values are scalars (number or quoted string), arrays (values in parentheses, separated by white space), tabular data.

### Object IDs
Object IDs are integers that uniquely identify each object in the file. They are used to link objects together in the file.

### Tables
Tables are a type of attribute that contain tabular data. They are used to store data that is not easily represented as a scalar or array.

