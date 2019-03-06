#pragma once

#include "LifeBrush.h"

// This stuff is registered in ShipEditor.cpp 
namespace FGraphVersion
{
	enum Type
	{
		Initial = 0,
		WithoutFGraphFace, // we removed graph faces (FGraphFace) from FGraph, also remove FGraphNode

		// -----------------------------------------
		// new versions can be added above this line
		VersionPlusOne,
		LatestVersion = VersionPlusOne - 1
	};

	const static FGuid GUID = FGuid( 0xD38E5AB2, 0x4782AFD7, 0xD96BB78E, 0xC8D4F792 );
};