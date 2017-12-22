// Fill out your copyright notice in the Description page of Project Settings.

using UnrealBuildTool;
using System.Collections.Generic;

public class MHPhysicsPrototype1Target : TargetRules
{
	public MHPhysicsPrototype1Target(TargetInfo Target) : base(Target)
	{
		Type = TargetType.Game;

		ExtraModuleNames.AddRange( new string[] { "MHPhysicsPrototype1" } );
	}
}
