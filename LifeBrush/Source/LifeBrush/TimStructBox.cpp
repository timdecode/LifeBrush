// Copyright 2018, Timothy Davison. All Rights Reserved.

#include "LifeBrush.h"

#include "Base64.h"

#include "TimStructBox.h"

void SkipWhitespace(const TCHAR*& Str)
{
	while (FChar::IsWhitespace(*Str))
	{
		Str++;
	}
}

bool SkipString(const TCHAR*& str, const FString& toSkip)
{
	const TCHAR * skip = *toSkip;

	const TCHAR * skipEnd = *toSkip + toSkip.Len();

	// find the start
	for (; *str != *skip; str++)
	{}

	// and skip
	for (; *str == *skip && skip < skipEnd; str++, skip++)
	{}

	if (skip != skipEnd)
		return false;

	return true;
}

void FTimStructBox::Destroy( UScriptStruct* ActualStruct )
{
	if(ActualStruct && structMemory)
	{
		ActualStruct->DestroyStruct( structMemory );
	}

	if(structMemory)
	{
		FMemory::Free( structMemory );
		structMemory = nullptr;
	}
}

void FTimStructBox::initStruct()
{
	if(!scriptStruct)
		return;

	Create( scriptStruct );
}

FTimStructBox::~FTimStructBox()
{
	ensure( scriptStruct || !structMemory );
	Destroy( scriptStruct );
}

FTimStructBox& FTimStructBox::operator=( const FTimStructBox& Other )
{
	if(this != &Other)
	{
		Destroy( scriptStruct );

		scriptStruct = Other.scriptStruct;

		if(Other.IsValid())
		{
			Create( scriptStruct );
			scriptStruct->CopyScriptStruct( structMemory, Other.structMemory );
		}
	}

	return *this;
}

void FTimStructBox::Create( UScriptStruct* ActualStruct )
{
	check( nullptr == structMemory );
	structMemory = (uint8*)FMemory::Malloc( ActualStruct->GetStructureSize() );
	ActualStruct->InitializeStruct( structMemory );
}

FTimStructBox::FTimStructBox( const FTimStructBox& Other )
{
	structMemory = nullptr;
	scriptStruct = Other.scriptStruct;

	if(Other.IsValid())
	{
		Create( scriptStruct );
		scriptStruct->CopyScriptStruct( structMemory, Other.structMemory );
	}
}

bool FTimStructBox::Serialize( FArchive& Ar )
{
	auto OldStruct = scriptStruct;
	Ar << scriptStruct;
	bool bValidBox = IsValid();
	Ar << bValidBox;

	if(Ar.IsLoading())
	{
		if(OldStruct != scriptStruct)
		{
			Destroy( OldStruct );
		}
		if(scriptStruct && !structMemory && bValidBox)
		{
			Create( scriptStruct );
		}
	}

	ensure( bValidBox || !IsValid() );
	if(IsValid() && bValidBox)
	{
		scriptStruct->SerializeItem( Ar, structMemory, nullptr );
	}

	return true;
}

// 

bool FTimStructBox::ExportTextItem(FString& ValueStr, FTimStructBox const& DefaultValue, UObject* Parent, int32 PortFlags, UObject* ExportRootScope) const
{
	ValueStr += TEXT("(");
	{
		UProperty * scriptStruct_property = FTimStructBox::StaticStruct()->FindPropertyByName(FName(TEXT("scriptStruct")));

		FString innerValue;
		if (scriptStruct_property->ExportText_InContainer(0, innerValue, this, &DefaultValue, Parent, PortFlags, ExportRootScope))
		{
			ValueStr += FString::Printf(TEXT("%s="), *scriptStruct_property->GetName());
			ValueStr += innerValue;
		}

		ValueStr += TEXT(",");

		ValueStr += TEXT("data=");

		scriptStruct->ExportText(ValueStr, structMemory, nullptr, Parent, PortFlags, ExportRootScope);

		//TArray<uint8> dataAsArray;

		//FMemoryWriter writer = FMemoryWriter(dataAsArray, true);

		//FObjectAndNameAsStringProxyArchive archive(writer, true);

		//// we're in a const function, yea, this could be dangerous, but we are writing not reading
		//FTimStructBox& self = *const_cast<FTimStructBox*>(this);
		//self.Serialize(archive);


		//ValueStr += FBase64::Encode(dataAsArray);

		//ValueStr += TEXT("'\"");

		//ValueStr += FBase64::Encode(structMemory, scriptStruct->GetStructureSize());
	}
	ValueStr += TEXT(")");

	return true;
}

//bool FTimStructBox::ExportTextItem(FString& ValueStr, FTimStructBox const& DefaultValue, UObject* Parent, int32 PortFlags, UObject* ExportRootScope) const
//{
//	ValueStr += TEXT("(");
//	{
//		UProperty * scriptStruct_property = FTimStructBox::StaticStruct()->FindPropertyByName(FName(TEXT("scriptStruct")));
//			
//		FString innerValue;
//		if (scriptStruct_property->ExportText_InContainer(0, innerValue, this, &DefaultValue, Parent, PortFlags, ExportRootScope))
//		{
//			ValueStr += FString::Printf(TEXT("%s="), *scriptStruct_property->GetName());
//			ValueStr += innerValue;
//		}
//
//		ValueStr += TEXT(",");
//
//		ValueStr += TEXT("data='\"");
//
//		TArray<uint8> dataAsArray;
//
//		FMemoryWriter writer = FMemoryWriter(dataAsArray, true);
//
//		FObjectAndNameAsStringProxyArchive archive(writer, true);
//
//		// we're in a const function, yea, this could be dangerous, but we are writing not reading
//		FTimStructBox& self = *const_cast<FTimStructBox*>(this);
//		self.Serialize(archive);
//
//
//		ValueStr += FBase64::Encode(dataAsArray);
//
//		ValueStr += TEXT("'\"");
//
//		//ValueStr += FBase64::Encode(structMemory, scriptStruct->GetStructureSize());
//	}
//	ValueStr += TEXT(")");
//
//	return true;
//}

/*

((scriptStruct=ScriptStruct'"/Script/RegionGrowing.FlexParticleObject"',data='"(inverseMass=0.125000)),(scriptStruct=ScriptStruct'"/Script/RegionGrowing.RandomWalkGraphObject"',data='"(maxVelocityOffset=1.000000,baseVelocity=1.000000)),(scriptStruct=ScriptStruct'"/Script/RegionGrowing.VelocityGraphObject"',data='"(linearVelocity=(X=0.000000,Y=0.000000,Z=0.000000),angularVelocity=(X=0.000000,Y=0.000000,Z=0.000000))),(scriptStruct=ScriptStruct'"/Script/RegionGrowing.ATPGraphObject"',data='"()),(scriptStruct=ScriptStruct'"/Script/RegionGrowing.GraphMesh"',data='"(staticMesh=StaticMesh'"/Game/2018_Cyberworlds/molecules/ims_a.ims_a"',material=MaterialInstanceConstant'"/Game/2018_Cyberworlds_Ext/materials/emissive_green_material.emissive_green_material"',visible=True)),(scriptStruct=ScriptStruct'"/Script/RegionGrowing.CachedElementGraphObject"',data='"(element=(radius=70.000000,generative=True,surfaceIndex=(faceIndex=-1,sectionIndex=-1),generationParameters=(),optimizationParameters=(),minAssignmentDistance=1.000000,freespaceRadius=1.000000,generationInnerRadius=1.000000,hackScale=1.000000))))
*/

/*

Begin Map
Begin Level
Begin Actor Class=/Script/Engine.Actor Name=a_simulation Archetype=/Script/Engine.Actor'/Script/Engine.Default__Actor'
Begin Object Class=/Script/RegionGrowing.FlexElements Name="FlexElements1"
Begin Object Class=/Script/RegionGrowing.GraphSimulationManager Name="GraphSimulationManager" Archetype=GraphSimulationManager'/Script/RegionGrowing.Default__FlexElements:GraphSimulationManager'
Begin Object Class=/Script/RegionGrowing.ATPSynthaseSimulation Name="ATPSynthaseSimulation_0"
End Object
Begin Object Class=/Script/RegionGrowing.TimelineSimulation Name="TimelineSimulation_0"
Begin Object Class=/Script/RegionGrowing.SEGraphTimeline Name="timeline" Archetype=SEGraphTimeline'/Script/RegionGrowing.Default__TimelineSimulation:timeline'
End Object
End Object
Begin Object Class=/Script/RegionGrowing.Visualization_AgentPathLines Name="Visualization_AgentPathLines_0"
End Object
End Object
End Object
Begin Object Class=/Script/Engine.SceneComponent Name="DefaultSceneRoot"
End Object
Begin Object Name="FlexElements1"
Begin Object Name="GraphSimulationManager"
Begin Object Name="ATPSynthaseSimulation_0"
atpTemplate(0)=(scriptStruct=ScriptStruct'"/Script/RegionGrowing.FlexParticleObject"',data=KQAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5GbGV4UGFydGljbGVPYmplY3QAAQAAAAYAAABncm91cAAMAAAASW50UHJvcGVydHkABAAAAAAAAAAAAAAAAAwAAABpbnZlcnNlTWFzcwAOAAAARmxvYXRQcm9wZXJ0eQAEAAAAAAAAAAAAAAA+CAAAAGNoYW5uZWwADAAAAEludFByb3BlcnR5AAQAAAAAAAAAAAAAAAAKAAAAbm9kZUluZGV4AAwAAABJbnRQcm9wZXJ0eQAEAAAAAAAAAAAAAAAABQAAAE5vbmUA)
atpTemplate(1)=(scriptStruct=ScriptStruct'"/Script/RegionGrowing.RandomWalkGraphObject"',data=LAAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5SYW5kb21XYWxrR3JhcGhPYmplY3QAAQAAAAkAAAB0aW1lTGVmdAAOAAAARmxvYXRQcm9wZXJ0eQAEAAAAAAAAAAAAAAAAEgAAAG1heFZlbG9jaXR5T2Zmc2V0AA4AAABGbG9hdFByb3BlcnR5AAQAAAAAAAAAAAAAgD8NAAAAYmFzZVZlbG9jaXR5AA4AAABGbG9hdFByb3BlcnR5AAQAAAAAAAAAAAAAgD8KAAAAbm9kZUluZGV4AAwAAABJbnRQcm9wZXJ0eQAEAAAAAAAAAAAAAAAABQAAAE5vbmUA)
atpTemplate(2)=(scriptStruct=ScriptStruct'"/Script/RegionGrowing.VelocityGraphObject"',data=KgAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5WZWxvY2l0eUdyYXBoT2JqZWN0AAEAAAAPAAAAbGluZWFyVmVsb2NpdHkADwAAAFN0cnVjdFByb3BlcnR5AAwAAAAAAAAABwAAAFZlY3RvcgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAAABhbmd1bGFyVmVsb2NpdHkADwAAAFN0cnVjdFByb3BlcnR5AAwAAAAAAAAABwAAAFZlY3RvcgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAoAAABub2RlSW5kZXgADAAAAEludFByb3BlcnR5AAQAAAAAAAAAAAAAAAAFAAAATm9uZQA=)
atpTemplate(3)=(scriptStruct=ScriptStruct'"/Script/RegionGrowing.ATPGraphObject"',data=JQAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5BVFBHcmFwaE9iamVjdAABAAAACgAAAG5vZGVJbmRleAAMAAAASW50UHJvcGVydHkABAAAAAAAAAAAAAAAAAUAAABOb25lAA==)
atpTemplate(4)=(scriptStruct=ScriptStruct'"/Script/RegionGrowing.GraphMesh"',data=IAAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5HcmFwaE1lc2gAAQAAAAsAAABzdGF0aWNNZXNoAA8AAABPYmplY3RQcm9wZXJ0eQAxAAAAAAAAAAAtAAAAL0dhbWUvMjAxOF9DeWJlcndvcmxkcy9tb2xlY3VsZXMvaW1zX2EuaW1zX2EACQAAAG1hdGVyaWFsAA8AAABPYmplY3RQcm9wZXJ0eQBZAAAAAAAAAABVAAAAL0dhbWUvMjAxOF9DeWJlcndvcmxkc19FeHQvbWF0ZXJpYWxzL2VtaXNzaXZlX2dyZWVuX21hdGVyaWFsLmVtaXNzaXZlX2dyZWVuX21hdGVyaWFsAAgAAAB2aXNpYmxlAA0AAABCb29sUHJvcGVydHkAAAAAAAAAAAABAAoAAABub2RlSW5kZXgADAAAAEludFByb3BlcnR5AAQAAAAAAAAAAAAAAAAFAAAATm9uZQA=)
atpTemplate(5)=(scriptStruct=ScriptStruct'"/Script/RegionGrowing.CachedElementGraphObject"',data=LwAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5DYWNoZWRFbGVtZW50R3JhcGhPYmplY3QAAQAAAAgAAABlbGVtZW50AA8AAABTdHJ1Y3RQcm9wZXJ0eQBtBAAAAAAAAAgAAABFbGVtZW50AAAAAAAAAAAAAAAAAAAAAAAABwAAAHJhZGl1cwAOAAAARmxvYXRQcm9wZXJ0eQAEAAAAAAAAAAAAAIxCCwAAAGVudGl0eVR5cGUADgAAAEludDE2UHJvcGVydHkAAgAAAAAAAAAAAAAFAAAAdHlwZQAOAAAASW50MTZQcm9wZXJ0eQACAAAAAAAAAAAAAAsAAABnZW5lcmF0aXZlAA0AAABCb29sUHJvcGVydHkAAAAAAAAAAAABAAwAAABlbnRpdHlJbmRleAAMAAAASW50UHJvcGVydHkABAAAAAAAAAAAAAAAAA0AAABzdXJmYWNlSW5kZXgADwAAAFN0cnVjdFByb3BlcnR5AGIAAAAAAAAADQAAAFN1cmZhY2VJbmRleAAAAAAAAAAAAAAAAAAAAAAAAAoAAABmYWNlSW5kZXgADAAAAEludFByb3BlcnR5AAQAAAAAAAAAAP////8NAAAAc2VjdGlvbkluZGV4AAwAAABJbnRQcm9wZXJ0eQAEAAAAAAAAAAD/////BQAAAE5vbmUAFQAAAGdlbmVyYXRpb25QYXJhbWV0ZXJzAA8AAABTdHJ1Y3RQcm9wZXJ0eQBdAAAAAAAAABgAAABOZWlnaGJvdXJob29kUGFyYW1ldGVycwAAAAAAAAAAAAAAAAAAAAAAAAcAAAByYWRpdXMADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAAAAAAkAAABrTmVhcmVzdAAMAAAASW50UHJvcGVydHkABAAAAAAAAAAAAAAAAAUAAABOb25lABcAAABvcHRpbWl6YXRpb25QYXJhbWV0ZXJzAA8AAABTdHJ1Y3RQcm9wZXJ0eQBdAAAAAAAAABgAAABOZWlnaGJvdXJob29kUGFyYW1ldGVycwAAAAAAAAAAAAAAAAAAAAAAAAcAAAByYWRpdXMADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAAAAAAkAAABrTmVhcmVzdAAMAAAASW50UHJvcGVydHkABAAAAAAAAAAAAAAAAAUAAABOb25lABYAAABtaW5Bc3NpZ25tZW50RGlzdGFuY2UADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAACAPxAAAABmcmVlc3BhY2VSYWRpdXMADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAACAPxYAAABnZW5lcmF0aW9uSW5uZXJSYWRpdXMADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAACAPwoAAABoYWNrU2NhbGUADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAACAPw0AAABncmFwaE9iamVjdHMADgAAAEFycmF5UHJvcGVydHkAUgAAAAAAAAAPAAAAU3RydWN0UHJvcGVydHkAAAAAAAANAAAAZ3JhcGhPYmplY3RzAA8AAABTdHJ1Y3RQcm9wZXJ0eQAAAAAAAAAAAA0AAABUaW1TdHJ1Y3RCb3gAAAAAAAAAAAAAAAAAAAAAAAAFAAAATm9uZQAKAAAAbm9kZUluZGV4AAwAAABJbnRQcm9wZXJ0eQAEAAAAAAAAAAAAAAAABQAAAE5vbmUA)
End Object
Begin Object Name="TimelineSimulation_0"
Begin Object Name="timeline"
End Object
defaultGlyph=(glyphMesh=StaticMesh'"/Engine/BasicShapes/Cube.Cube"',glyphMaterial=MaterialInstanceConstant'"/Game/SIGGRAPH_Asia2017/basic_materials/dark_dark_grey.dark_dark_grey"',eventClass=Class'"/Script/RegionGrowing.SEGraphEvent"')
glyphPrototypes(0)=(glyphMesh=StaticMesh'"/Engine/BasicShapes/Cube.Cube"',glyphMaterial=MaterialInstanceConstant'"/Game/2018_Cyberworlds_Ext/materials/emissive_red_material.emissive_red_material"',eventClass=Class'"/Script/RegionGrowing.Event_ProtonPumped"')
glyphPrototypes(1)=(glyphMesh=StaticMesh'"/Engine/BasicShapes/Cube.Cube"',glyphMaterial=MaterialInstanceConstant'"/Game/2018_Cyberworlds_Ext/materials/emissive_green_material.emissive_green_material"',eventClass=Class'"/Script/RegionGrowing.Event_SpawnATP"')
glyphPrototypes(2)=(glyphMesh=StaticMesh'"/Engine/BasicShapes/Cube.Cube"',glyphMaterial=Material'"/Game/2018_Cyberworlds_Ext/materials/lineCompatible_emissiveMaterial.lineCompatible_emissiveMaterial"',eventClass=Class'"/Script/RegionGrowing.Event_SpinATPSynthase"')
timeline="timeline"
End Object
Begin Object Name="Visualization_AgentPathLines_0"
defaultLineMaterial=Material'"/Game/DebugRed.DebugRed"'
End Object
simulations(0)=ATPSynthaseSimulation'"ATPSynthaseSimulation_0"'
simulations(1)=TimelineSimulation'"TimelineSimulation_0"'
simulations(2)=Visualization_AgentPathLines'"Visualization_AgentPathLines_0"'
End Object
flexParams=(gravity=(X=0.000000,Y=0.000000,Z=0.000000),radius=7.000000,solidRestDistance=1.800000,fluidRestDistance=2.000000,restitution=0.200000,dissipation=0.010000,cohesion=0.010000,viscosity=0.010000,buoyancy=0.000000,collisionDistance=1.000000,particleCollisionMargin=1.000000,shapeCollisionMargin=1.000000,planes[0]=(X=1.000000,Y=0.000000,Z=0.000000,W=80.738922),planes[1]=(X=-1.000000,Y=0.000000,Z=0.000000,W=3.011204),planes[2]=(X=0.000000,Y=1.000000,Z=0.000000,W=-271.246979),planes[3]=(X=0.000000,Y=-1.000000,Z=0.000000,W=424.488373),planes[4]=(X=0.000000,Y=0.000000,Z=1.000000,W=-37.346622),planes[5]=(X=0.000000,Y=0.000000,Z=-1.000000,W=149.079468),numPlanes=6,relaxationMode=0)
simulationManager="GraphSimulationManager"
CreationMethod=Instance
End Object
Begin Object Name="DefaultSceneRoot"
RelativeLocation=(X=-8.306232,Y=-232.057602,Z=36.511074)
bVisualizeComponent=True
CreationMethod=Instance
End Object
bNetStartup=True
RootComponent=SceneComponent'"DefaultSceneRoot"'
ActorLabel="a_simulation"
bIsEditorPreviewActor=True
InstanceComponents(0)=SceneComponent'"DefaultSceneRoot"'
InstanceComponents(1)=FlexElements'"FlexElements1"'
End Actor
End Level
Begin Surface
End Surface
End Map


*/


/*

((scriptStruct=ScriptStruct'"/Script/RegionGrowing.FlexParticleObject"',data='"KQAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5GbGV4UGFydGljbGVPYmplY3QAAQAAAAYAAABncm91cAAMAAAASW50UHJvcGVydHkABAAAAAAAAAAAAAAAAAwAAABpbnZlcnNlTWFzcwAOAAAARmxvYXRQcm9wZXJ0eQAEAAAAAAAAAAAAAAA+CAAAAGNoYW5uZWwADAAAAEludFByb3BlcnR5AAQAAAAAAAAAAAAAAAAKAAAAbm9kZUluZGV4AAwAAABJbnRQcm9wZXJ0eQAEAAAAAAAAAAAAAAAABQAAAE5vbmUA'"),(scriptStruct=ScriptStruct'"/Script/RegionGrowing.RandomWalkGraphObject"',data='"LAAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5SYW5kb21XYWxrR3JhcGhPYmplY3QAAQAAAAkAAAB0aW1lTGVmdAAOAAAARmxvYXRQcm9wZXJ0eQAEAAAAAAAAAAAAAAAAEgAAAG1heFZlbG9jaXR5T2Zmc2V0AA4AAABGbG9hdFByb3BlcnR5AAQAAAAAAAAAAAAAgD8NAAAAYmFzZVZlbG9jaXR5AA4AAABGbG9hdFByb3BlcnR5AAQAAAAAAAAAAAAAgD8KAAAAbm9kZUluZGV4AAwAAABJbnRQcm9wZXJ0eQAEAAAAAAAAAAAAAAAABQAAAE5vbmUA'"),(scriptStruct=ScriptStruct'"/Script/RegionGrowing.VelocityGraphObject"',data='"KgAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5WZWxvY2l0eUdyYXBoT2JqZWN0AAEAAAAPAAAAbGluZWFyVmVsb2NpdHkADwAAAFN0cnVjdFByb3BlcnR5AAwAAAAAAAAABwAAAFZlY3RvcgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAAABhbmd1bGFyVmVsb2NpdHkADwAAAFN0cnVjdFByb3BlcnR5AAwAAAAAAAAABwAAAFZlY3RvcgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAoAAABub2RlSW5kZXgADAAAAEludFByb3BlcnR5AAQAAAAAAAAAAAAAAAAFAAAATm9uZQA='"),(scriptStruct=ScriptStruct'"/Script/RegionGrowing.ATPGraphObject"',data='"JQAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5BVFBHcmFwaE9iamVjdAABAAAACgAAAG5vZGVJbmRleAAMAAAASW50UHJvcGVydHkABAAAAAAAAAAAAAAAAAUAAABOb25lAA=='"),(scriptStruct=ScriptStruct'"/Script/RegionGrowing.GraphMesh"',data='"IAAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5HcmFwaE1lc2gAAQAAAAsAAABzdGF0aWNNZXNoAA8AAABPYmplY3RQcm9wZXJ0eQAxAAAAAAAAAAAtAAAAL0dhbWUvMjAxOF9DeWJlcndvcmxkcy9tb2xlY3VsZXMvaW1zX2EuaW1zX2EACQAAAG1hdGVyaWFsAA8AAABPYmplY3RQcm9wZXJ0eQBZAAAAAAAAAABVAAAAL0dhbWUvMjAxOF9DeWJlcndvcmxkc19FeHQvbWF0ZXJpYWxzL2VtaXNzaXZlX2dyZWVuX21hdGVyaWFsLmVtaXNzaXZlX2dyZWVuX21hdGVyaWFsAAgAAAB2aXNpYmxlAA0AAABCb29sUHJvcGVydHkAAAAAAAAAAAABAAoAAABub2RlSW5kZXgADAAAAEludFByb3BlcnR5AAQAAAAAAAAAAAAAAAAFAAAATm9uZQA='"),(scriptStruct=ScriptStruct'"/Script/RegionGrowing.CachedElementGraphObject"',data='"LwAAAC9TY3JpcHQvUmVnaW9uR3Jvd2luZy5DYWNoZWRFbGVtZW50R3JhcGhPYmplY3QAAQAAAAgAAABlbGVtZW50AA8AAABTdHJ1Y3RQcm9wZXJ0eQBtBAAAAAAAAAgAAABFbGVtZW50AAAAAAAAAAAAAAAAAAAAAAAABwAAAHJhZGl1cwAOAAAARmxvYXRQcm9wZXJ0eQAEAAAAAAAAAAAAAIxCCwAAAGVudGl0eVR5cGUADgAAAEludDE2UHJvcGVydHkAAgAAAAAAAAAAAAAFAAAAdHlwZQAOAAAASW50MTZQcm9wZXJ0eQACAAAAAAAAAAAAAAsAAABnZW5lcmF0aXZlAA0AAABCb29sUHJvcGVydHkAAAAAAAAAAAABAAwAAABlbnRpdHlJbmRleAAMAAAASW50UHJvcGVydHkABAAAAAAAAAAAAAAAAA0AAABzdXJmYWNlSW5kZXgADwAAAFN0cnVjdFByb3BlcnR5AGIAAAAAAAAADQAAAFN1cmZhY2VJbmRleAAAAAAAAAAAAAAAAAAAAAAAAAoAAABmYWNlSW5kZXgADAAAAEludFByb3BlcnR5AAQAAAAAAAAAAP////8NAAAAc2VjdGlvbkluZGV4AAwAAABJbnRQcm9wZXJ0eQAEAAAAAAAAAAD/////BQAAAE5vbmUAFQAAAGdlbmVyYXRpb25QYXJhbWV0ZXJzAA8AAABTdHJ1Y3RQcm9wZXJ0eQBdAAAAAAAAABgAAABOZWlnaGJvdXJob29kUGFyYW1ldGVycwAAAAAAAAAAAAAAAAAAAAAAAAcAAAByYWRpdXMADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAAAAAAkAAABrTmVhcmVzdAAMAAAASW50UHJvcGVydHkABAAAAAAAAAAAAAAAAAUAAABOb25lABcAAABvcHRpbWl6YXRpb25QYXJhbWV0ZXJzAA8AAABTdHJ1Y3RQcm9wZXJ0eQBdAAAAAAAAABgAAABOZWlnaGJvdXJob29kUGFyYW1ldGVycwAAAAAAAAAAAAAAAAAAAAAAAAcAAAByYWRpdXMADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAAAAAAkAAABrTmVhcmVzdAAMAAAASW50UHJvcGVydHkABAAAAAAAAAAAAAAAAAUAAABOb25lABYAAABtaW5Bc3NpZ25tZW50RGlzdGFuY2UADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAACAPxAAAABmcmVlc3BhY2VSYWRpdXMADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAACAPxYAAABnZW5lcmF0aW9uSW5uZXJSYWRpdXMADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAACAPwoAAABoYWNrU2NhbGUADgAAAEZsb2F0UHJvcGVydHkABAAAAAAAAAAAAACAPw0AAABncmFwaE9iamVjdHMADgAAAEFycmF5UHJvcGVydHkAUgAAAAAAAAAPAAAAU3RydWN0UHJvcGVydHkAAAAAAAANAAAAZ3JhcGhPYmplY3RzAA8AAABTdHJ1Y3RQcm9wZXJ0eQAAAAAAAAAAAA0AAABUaW1TdHJ1Y3RCb3gAAAAAAAAAAAAAAAAAAAAAAAAFAAAATm9uZQAKAAAAbm9kZUluZGV4AAwAAABJbnRQcm9wZXJ0eQAEAAAAAAAAAAAAAAAABQAAAE5vbmUA'"))

*/

bool FTimStructBox::ImportTextItem(const TCHAR*& InBuffer, int32 PortFlags, UObject* Parent, FOutputDevice* ErrorText)
{
	TArray<FDefinedProperty> DefinedProperties;

	FString StructName = FTimStructBox::StaticStruct()->GetName();

	TArray<uint8> data;

	// this keeps track of the number of errors we've logged, so that we can add new lines when logging more than one error
	int32 ErrorCount = 0;
	const TCHAR* Buffer = InBuffer;
	if (*Buffer++ == TCHAR('('))
	{
		// parse and import the value
		Buffer = UProperty::ImportSingleProperty(Buffer, this, FTimStructBox::StaticStruct(), Parent, PortFlags | PPF_Delimited, ErrorText, DefinedProperties);

		// Buffer is sitting at the end of the property name
		// skip any remaining text before the next property value
		{
			SkipWhitespace(Buffer);
			int32 SubCount = 0;
			while (*Buffer && *Buffer != TCHAR('\r') && *Buffer != TCHAR('\n') &&
				(SubCount > 0 || *Buffer != TCHAR(')')) && (SubCount > 0 || *Buffer != TCHAR(',')))
			{
				SkipWhitespace(Buffer);
				if (*Buffer == TCHAR('\"'))
				{
					do
					{
						Buffer++;
					} while (*Buffer && *Buffer != TCHAR('\"') && *Buffer != TCHAR('\n') && *Buffer != TCHAR('\r'));

					if (*Buffer != TCHAR('\"'))
					{
						ErrorText->Logf(TEXT("%sImportText (%s): Bad quoted string at: %s"), ErrorCount++ > 0 ? LINE_TERMINATOR : TEXT(""), *StructName, Buffer);
						return nullptr;
					}
				}
				else if (*Buffer == TCHAR('('))
				{
					SubCount++;
				}
				else if (*Buffer == TCHAR(')'))
				{
					SubCount--;
					if (SubCount < 0)
					{
						ErrorText->Logf(TEXT("%sImportText (%s): Too many closing parenthesis in: %s"), ErrorCount++ > 0 ? LINE_TERMINATOR : TEXT(""), *StructName, InBuffer);
						return nullptr;
					}
				}
				Buffer++;
			}
			if (SubCount > 0)
			{
				ErrorText->Logf(TEXT("%sImportText(%s): Not enough closing parenthesis in: %s"), ErrorCount++ > 0 ? LINE_TERMINATOR : TEXT(""), *StructName, InBuffer);
				return nullptr;
			}
		}



		// Skip comma
		if (*Buffer == TCHAR(','))
		{
			// Skip comma
			Buffer++;
		}

		if (scriptStruct)
		{
			Destroy(scriptStruct);
		}

		if (scriptStruct && !structMemory)
		{
			Create(scriptStruct);
		}

		const FString dataHeader(TEXT("data="));
		SkipString(Buffer, dataHeader);

		Buffer = scriptStruct->ImportText(Buffer, structMemory, Parent, PortFlags, ErrorText, scriptStruct->GetName());


		
		// sit at the end
		Buffer = Buffer + 1;
	}

	if (DefinedProperties.Num() != 1)
		return false;


	// all good, set the InBuffer
	InBuffer = Buffer;


	return true;
}

bool FTimStructBox::Identical( const FTimStructBox* Other, uint32 PortFlags ) const
{
	if(!Other)
	{
		return false;
	}

	if(scriptStruct != Other->scriptStruct)
	{
		return false;
	}

	if(!scriptStruct)
	{
		return true;
	}

	if(!structMemory && !Other->structMemory)
	{
		return true;
	}

	return scriptStruct->CompareScriptStruct( structMemory, Other->structMemory, PortFlags );
}

void FTimStructBox::AddStructReferencedObjects( class FReferenceCollector& Collector ) const
{
	Collector.AddReferencedObject( const_cast<FTimStructBox*>(this)->scriptStruct );
	if(scriptStruct && structMemory && scriptStruct->RefLink)
	{
		scriptStruct->SerializeBin( Collector.GetVerySlowReferenceCollectorArchive(), structMemory );
	}
}
