# Style Guide
## General
- Line length should not exceed 100 characters
- Function bodies should not exceed 40 lines of code (excluding comments and whitespace)

## License
- A license should appear at the top of every file

## Layout
### Defines
- Use the `#ifndef` - `#define` paradigm instead of `#pragma once`

- `#define` preprocessor macros are defined before `#include`s
- Use preprocessor macros for compilation flags and constants
- Avoid using preprocessor macros as functions - where possible, use a real (potentially inlined) function instead

### Includes
- In `.cc` files, the corresponding `.h` file should always be included first
- Prefer quotations instead of angle brackets when including project files
- System files (`#include <...>`) are included before other project files (`#include "..."`)
- Alphabetical ordering within each group of `#include`s

### Class definitions
- All required forward declarations are declared before the class definition

- Friends declared at the end of the class
- Inline functions defined after the class definition - NOT within the class

- Types (internal structs, typedefs, etc.) are declared before variables, which are declared before functions
- `private` types are declared before `protected` types, which are declared before `public` types
- `private` variables are declared before `protected` variables, which are declared before `public` variables
- Constructors are declared before all other functions, regardless of visibility
- `public` functions are declared before `protected` functions, which are declared before `private` functions
- `static` members are declared before non-`static` members

- Within each `public`/`protected`/`private` block, methods should be declared in the following order: accessors, mutators, statistics, heuristics, event handlers, helper functions

### Function signatures
- Prefer to return values instead of modifying parameters
- Modified parameters occur earlier in the parameter list

### Function bodies
- Within a function, code should not be nested more than three times
    - Adjust the control flow to exit early and extract functionality to helper functions

## Naming
- PascalCase for type definitions
- camelCase for variables
- SNAKE_CAPS for constants
- avoid snake_case

- File names should not exceed 64 characters
- Prefer descriptive names over concise names
- Don't abbreviate variable names

## Brace style
- [1TBS](https://en.wikipedia.org/wiki/Indentation_style#Variant:_1TBS_(OTBS))
- Closing braces occur on their own line

## Spacing and Indentation
- 1 tab = 4 spaces
- Code within curly braces is indented one level
- No space between function name and round bracket

### Defines
- Where possible, `#if`s should be backward-indented once

### Keywords
- `public`, `protected`, and `private` should all be indented at the same level as the initial class/struct definition
- `case` and `default` should be indented at the same level as `switch`

## Comments
- Every member type has doxygen-style comments when the type is defined
- Every member variable has doxygen-style comments when the variable is declared
- Every member function has doxygen-style comments when the signature is declared
- Every class has doxygen-style comments when the class is defined

- Functions performing multi-step processes should have a comment before each step, describing the step