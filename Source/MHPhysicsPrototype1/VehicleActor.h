#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Pawn.h"
#include "MHPhysics/MHStaticMeshComponent.h"
#include "VehicleActor.generated.h"

USTRUCT()
struct FVehicleControl
{
	GENERATED_BODY()

public:
	FVehicleControl();

	void Tick(float DeltaSeconds);

	void SetThrottle(float InThrottle);
	void SetBrake(float InBrake);
	void SetSteer(float InSteer);

	float GetThrottle() const { return Throttle; }
	float GetBrake() const { return Brake; }
	float GetSteer() const { return Steer; }
	int32 GetGear() const { return Gear; }

	void OnThrottle(float InThrottle);
	void OnBrake(float InBrake);
	void OnSteer(float InSteer);
	void OnSteerLeftPressed();
	void OnSteerLeftReleased();
	void OnSteerRightPressed();
	void OnSteerRightReleased();

private:
	void TickButtonSteering(float DeltaSeconds);

	UPROPERTY(EditAnywhere, Category = "Vehicle Input")
	float Throttle;

	UPROPERTY(EditAnywhere, Category = "Vehicle Input")
	float Brake;

	UPROPERTY(EditAnywhere, Category = "Vehicle Input")
	float Steer;

	UPROPERTY(EditAnywhere, Category = "Vehicle Input")
	int32 Gear;

	bool bSteerByButton;
	float SteerAxisInput;
	bool bSteerLeft;
	bool bSteerRight;
};


UCLASS()
class MHPHYSICSPROTOTYPE1_API AVehicleActor : public APawn
{
	GENERATED_BODY()

public:
	// Sets default values for this pawn's properties
	AVehicleActor();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	// Called to bind functionality to input
	virtual void SetupPlayerInputComponent(class UInputComponent* PlayerInputComponent) override;

private:
	UPROPERTY()
	FVehicleControl Control;

	UPROPERTY()
	UMHStaticMeshComponent* MeshComponent;
	
};
