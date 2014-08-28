#define SECTION_LIBRARY 0
#define SECTION_CONJUGATES 1
#define SECTION_RULES 2
#define SECTION_INITIAL 3
#define SECTION_CONTROL 4
#define SECTION_BOUNDARY 5
#define SECTION_GEOMETRY 6
#define SECTION_END 7
#define SECTION_CONFIG 8

extern ChemistryComputation parseConfigFile(ifstream& file);
