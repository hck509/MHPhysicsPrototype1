// Fill out your copyright notice in the Description page of Project Settings.

#include "MHPhysicsPrototype1GameModeBase.h"
#include "VehicleActor.h"
#include "EngineUtils.h"

AMHPhysicsPrototype1GameModeBase::AMHPhysicsPrototype1GameModeBase()
{
	PrimaryActorTick.bCanEverTick = true;
	PrimaryActorTick.bStartWithTickEnabled = true;
}

void AMHPhysicsPrototype1GameModeBase::StartPlay()
{
	Super::StartPlay();

	MHPhysics.GenerateFromStaticMesheActors(GetWorld());

	for (TActorIterator<AVehicleActor> It(GetWorld()); It; ++It)
	{
		AVehicleActor* Vehicle = *It;
		GetWorld()->GetFirstPlayerController()->Possess(Vehicle);
		break;
	}
}

void AMHPhysicsPrototype1GameModeBase::Tick(float DeltaSeconds)
{
	Super::Tick(DeltaSeconds);

	MHPhysics.Tick(DeltaSeconds);
	MHPhysics.DebugDraw(GetWorld());
}
