// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Components/MeshComponent.h"
#include "MHPhysics/MHPHysics.h"
#include "MHStaticMeshComponent.generated.h"


UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class MHPHYSICSPROTOTYPE1_API UMHStaticMeshComponent : public UMeshComponent
{
	GENERATED_BODY()

public:	
	// Sets default values for this component's properties
	UMHStaticMeshComponent();

protected:
	// Called when the game starts
	virtual void BeginPlay() override;
	virtual void OnRegister() override;

public:	
	// Called every frame
	virtual void TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) override;

	void InitializeCustomMesh();
		
private:
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = StaticMesh, meta = (AllowPrivateAccess = "true"))
	class UStaticMesh* StaticMesh;

	UPROPERTY(BlueprintReadOnly, Category = Mesh, meta = (AllowPrivateAccess = "true"))
	class UCustomMeshComponent* CustomMeshComponent;

	UPROPERTY(EditInstanceOnly, BlueprintReadOnly, Category = MHPHysics, meta = (AllowPrivateAccess = "true"))
	float MassInKg;

	UPROPERTY(EditInstanceOnly, BlueprintReadOnly, Category = MHPHysics, meta = (AllowPrivateAccess = "true"))
	float SpringK;

	UPROPERTY(EditInstanceOnly, BlueprintReadOnly, Category = MHPHysics, meta = (AllowPrivateAccess = "true"))
	float SpringD;

	FMHMeshInfo MHMeshInfo;
};
