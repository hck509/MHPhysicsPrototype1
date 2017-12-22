// Fill out your copyright notice in the Description page of Project Settings.

#include "MHPhysicsPrototype1GameModeBase.h"


AMHPhysicsPrototype1GameModeBase::AMHPhysicsPrototype1GameModeBase()
{
	PrimaryActorTick.bCanEverTick = true;
	PrimaryActorTick.bStartWithTickEnabled = true;
}

void AMHPhysicsPrototype1GameModeBase::StartPlay()
{
	Super::StartPlay();

	MHPhysics.GenerateFromStaticMesheActors(GetWorld());
}

void AMHPhysicsPrototype1GameModeBase::Tick(float DeltaSeconds)
{
	Super::Tick(DeltaSeconds);

	MHPhysics.Tick(DeltaSeconds);
	MHPhysics.DebugDraw(GetWorld());
}
