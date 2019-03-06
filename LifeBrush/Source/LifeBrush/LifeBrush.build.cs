// Fill out your copyright notice in the Description page of Project Settings.

using UnrealBuildTool;

using System;
using System.IO;
using System.Collections.Generic;

public class LifeBrush : ModuleRules
{
    public LifeBrush(ReadOnlyTargetRules Target) : base(Target)
	{
		PublicDependencyModuleNames.AddRange(new string[] {
            "Core",
            "CoreUObject",
            "Engine",
            "InputCore",
            "RenderCore",
            "UMG",
            "Slate",
            "RHI",
            "ShaderCore",
            "RuntimeMeshComponent"
        });
        PrivateDependencyModuleNames.AddRange(new string[] {
            "Slate",
            "SlateCore",
            "InputCore",
            "HeadMountedDisplay"
        });

        if (Target.bBuildEditor)
        {
            PublicDependencyModuleNames.AddRange(new string[] {"UnrealEd", "LevelEditor"});

            PrivateDependencyModuleNames.AddRange(
                new string[] {
                    "EditorStyle",
                    "UnrealEd",
                    "BspMode",
                    "LevelEditor",
                    "PropertyEditor",
                    "BspMode"
                }
            );
        }
        
        _addTcods();
        _addFlex();
        _addEigen();

        
        bEnableShadowVariableWarnings = false;

        MinFilesUsingPrecompiledHeaderOverride = 1;
        bFasterWithoutUnity = true;
	}

    private void _addEigen()
    {
        PublicIncludePaths.Add(Path.Combine(ModuleDirectory, "eigenlib"));
    }

    private void _addFlex()
    {
        // add the .libs
        string libraryPath = Path.Combine(ModuleDirectory, "Flex", "lib");

        string platformPath = "";

        string[] libraries = {};

        if (Target.Platform == UnrealTargetPlatform.Win64)
        {
            platformPath = "Win64";

            libraries = new string[] {
                "NvFlexReleaseCUDA_x64.lib",
                "NvFlexReleaseD3D_x64.lib",
                "NvFlexExtReleaseCUDA_x64.lib",
                "NvFlexExtReleaseD3D_x64.lib",
                "NvFlexDeviceRelease_x64.lib"
            };
        }

        libraryPath = Path.Combine(libraryPath, platformPath);

        List<string> libraryPaths = new List<string>();
        
        foreach( string lib in libraries )
        {
           libraryPaths.Add(Path.Combine(libraryPath, lib));
        }

        PublicAdditionalLibraries.AddRange(libraryPaths.ToArray());

        // add the include
        PublicIncludePaths.Add(Path.Combine(ModuleDirectory, "Flex", "include"));
    }

    private void _addTcods()
    {
        string libraryPath = Path.Combine(ModuleDirectory, "Algorithm", "tcods", "lib");

        string platformPath = "";

        string[] libraries = {};

        if (Target.Platform == UnrealTargetPlatform.Win64)
        {
            platformPath = "Win64";
            libraries = new string[] {
                "suitesparseconfig.lib",
                "libamd.lib",
                "libcamd.lib",
                "libcolamd.lib",
                "libccolamd.lib",
                "libcholmod.lib",
                "libspqr.lib",
                "metis.lib",
                "liblapack.lib",
                "libblas.lib",
                "libumfpack.lib"
            };
        }
        else if (Target.Platform == UnrealTargetPlatform.Mac)
        {
            platformPath = "Mac";

            libraries = new string[] {
                "libsuitesparseconfig.a",
                "libamd.a",
                "libcamd.a",
                "libcolamd.a",
                "libccolamd.a",
                "libcholmod.a",
                "libspqr.a",
                "libmetis.a",
                "libumfpack.a"
            };
        }

        libraryPath = Path.Combine(libraryPath, platformPath);

        List<string> libraryPaths = new List<string>();
        
        foreach( string lib in libraries )
        {
           libraryPaths.Add(Path.Combine(libraryPath, lib));
        }

        if (Target.Platform == UnrealTargetPlatform.Mac)
        {
            libraryPaths.Add("m");
            PublicAdditionalFrameworks.Add(new UEBuildFramework("Accelerate"));
        }

        PublicAdditionalLibraries.AddRange(libraryPaths.ToArray());
        PublicIncludePaths.Add(Path.Combine(ModuleDirectory, "Algorithm", "tcods", "include"));
    }
}
