// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"

class FMHPhysicsPrototype1GameModule : public FDefaultGameModuleImpl
{
public:
	FMHPhysicsPrototype1GameModule();

	virtual void StartupModule() override;
	virtual void ShutdownModule() override;
};
