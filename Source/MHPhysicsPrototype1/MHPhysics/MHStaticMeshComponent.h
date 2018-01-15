// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Components/MeshComponent.h"
#include "MHPhysics/MHPHysics.h"

#if WITH_EDITORONLY_DATA
#include "IDetailCustomization.h"
#endif

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

	const FMHMeshInfo& GetMHMeshInfo() const { return MHMeshInfo; }

#if WITH_EDITORONLY_DATA
	void ImportFBX();
	void Reimport();
#endif
		
private:
	void InitializeFromStaticMesh();
	void InitializeFromChunk();

	void UpdateCustomMeshFromChunk();
	void UpdateCustomMeshFromMHPhysics();

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

	UPROPERTY()
	FMHChunk MHChunk;

	FMHMeshInfo MHMeshInfo;

#if WITH_EDITORONLY_DATA
	UPROPERTY()
	FString SourceFilePath;
#endif
};

#if WITH_EDITORONLY_DATA

class FMHStaticMeshComponentDetails : public IDetailCustomization
{
public:
	/** Makes a new instance of this detail layout class for a specific detail view requesting it */
	static TSharedRef<IDetailCustomization> MakeInstance();

private:
	/** IDetailCustomization interface */
	virtual void CustomizeDetails(IDetailLayoutBuilder& DetailLayout) override;

	FReply ImportFBX();
	FReply Reimport();

private:
	TArray<TWeakObjectPtr<UMHStaticMeshComponent>> Components;
};

#endif