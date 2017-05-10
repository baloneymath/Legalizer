ADD_LIBRARY (Legalizer
Legal.cpp
Abacus.cpp
)
INCLUDE_DIRECTORIES (
        ${PROJECT_SOURCE_DIR}/Parser
        ${PROJECT_SOURCE_DIR}/PlaceCommon
        ${PROJECT_SOURCE_DIR}/Placement
	${PROJECT_SOURCE_DIR}/Gnuplot
	${PROJECT_SOURCE_DIR}/
)
