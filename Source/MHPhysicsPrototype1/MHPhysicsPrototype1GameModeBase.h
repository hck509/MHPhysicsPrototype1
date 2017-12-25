// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/GameModeBase.h"
#include "MHPhysics/MHPHysics.h"
#include "MHPhysicsPrototype1GameModeBase.generated.h"

/**
 * 
 */
UCLASS()
class MHPHYSICSPROTOTYPE1_API AMHPhysicsPrototype1GameModeBase : public AGameModeBase
{
	GENERATED_BODY()
public:

	AMHPhysicsPrototype1GameModeBase();
	
	virtual void StartPlay() override;
	virtual void Tick(float DeltaSeconds) override;

	FMHPhysics& GetMHPhyscis() { return MHPhysics; }

private:
	FMHPhysics MHPhysics;
};
