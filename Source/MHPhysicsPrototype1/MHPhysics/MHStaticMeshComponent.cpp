// Fill out your copyright notice in the Description page of Project Settings.

#include "MHStaticMeshComponent.h"
#include "CustomMeshComponent.h"
#include "MHPhysicsPrototype1GameModeBase.h"
#include "Engine/World.h"

// Sets default values for this component's properties
UMHStaticMeshComponent::UMHStaticMeshComponent()
{
	// Set this component to be initialized when the game starts, and to be ticked every frame.  You can turn these features
	// off to improve performance if you don't need them.
	PrimaryComponentTick.bCanEverTick = true;

	CustomMeshComponent = CreateDefaultSubobject<UCustomMeshComponent>(TEXT("CustomMeshComponent"));
	CustomMeshComponent->SetupAttachment(this);

	MassInKg = 100.0f;
	SpringK = 1000.0f;
	SpringD = 10.0f;
}

// Called when the game starts
void UMHStaticMeshComponent::BeginPlay()
{
	Super::BeginPlay();
}

void UMHStaticMeshComponent::OnRegister()
{
	Super::OnRegister();

	InitializeCustomMesh();
}

// Called every frame
void UMHStaticMeshComponent::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

}

void UMHStaticMeshComponent::InitializeCustomMesh()
{
	if (StaticMesh)
	{
		if (ensure(StaticMesh->RenderData.IsValid()) && ensure(StaticMesh->RenderData->LODResources.Num() > 0))
		{
			TArray<FCustomMeshTriangle> Triangles;

			const FTransform ComponentTransform = GetComponentTransform();

			const FStaticMeshLODResources& Resource = StaticMesh->RenderData->LODResources[0];

			FIndexArrayView IndexArrayView = Resource.IndexBuffer.GetArrayView();

			ensure(IndexArrayView.Num() % 3 == 0);

			for (int32 Index = 0; Index + 2 < IndexArrayView.Num(); Index += 3)
			{
				Triangles.Add(FCustomMeshTriangle({
					Resource.PositionVertexBuffer.VertexPosition(IndexArrayView[Index]),
					Resource.PositionVertexBuffer.VertexPosition(IndexArrayView[Index + 1]),
					Resource.PositionVertexBuffer.VertexPosition(IndexArrayView[Index + 2]) }));
			}

			CustomMeshComponent->SetCustomMeshTriangles(Triangles);

			AMHPhysicsPrototype1GameModeBase* GameMode = Cast<AMHPhysicsPrototype1GameModeBase>(GetWorld()->GetAuthGameMode());

			if (GameMode)
			{
				MHMeshInfo = GameMode->GetMHPhyscis().GenerateFromStaticMesh(*StaticMesh, GetComponentTransform(), MassInKg, SpringK, SpringD);
			}
		}
	}
}

