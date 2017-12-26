// Fill out your copyright notice in the Description page of Project Settings.

#include "MHPhysicsPrototype1.h"
#include "Modules/ModuleManager.h"
#include "MHPhysics/MHStaticMeshComponent.h"
#include "PropertyEditorModule.h"

IMPLEMENT_PRIMARY_GAME_MODULE(FMHPhysicsPrototype1GameModule, MHPhysicsPrototype1, "MHPhysicsPrototype1");

FMHPhysicsPrototype1GameModule::FMHPhysicsPrototype1GameModule()
{

}

void FMHPhysicsPrototype1GameModule::StartupModule()
{
	FPropertyEditorModule& PropertyModule = FModuleManager::LoadModuleChecked<FPropertyEditorModule>("PropertyEditor");
	{
		PropertyModule.RegisterCustomClassLayout(UMHStaticMeshComponent::StaticClass()->GetFName(), FOnGetDetailCustomizationInstance::CreateStatic(&FMHStaticMeshComponentDetails::MakeInstance));
	}
}

void FMHPhysicsPrototype1GameModule::ShutdownModule()
{
	FPropertyEditorModule& PropertyModule = FModuleManager::LoadModuleChecked<FPropertyEditorModule>("PropertyEditor");

	PropertyModule.UnregisterCustomClassLayout(UMHStaticMeshComponent::StaticClass()->GetFName());
}

