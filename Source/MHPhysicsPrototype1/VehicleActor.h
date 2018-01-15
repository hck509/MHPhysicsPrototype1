#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Pawn.h"
#include "MHPhysics/MHStaticMeshComponent.h"
#include "VehicleActor.generated.h"

UCLASS()
class UVehicleControl : public UObject
{
	GENERATED_BODY()

public:
	UVehicleControl();

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

	virtual void PossessedBy(AController* NewController) override;
	virtual void SetupPlayerInputComponent(class UInputComponent* PlayerInputComponent) override;

private:
	UPROPERTY()
	UVehicleControl* Control;

	UPROPERTY(EditAnywhere)
	UMHStaticMeshComponent* MeshComponent;
	
	/** Spring arm that will offset the camera */
	UPROPERTY(Category = "Vehicle Camera", VisibleDefaultsOnly, BlueprintReadOnly, meta = (AllowPrivateAccess = "true"))
	class USpringArmComponent* CameraSpringArm;

	/** Camera component that will be our viewpoint */
	UPROPERTY(Category = "Vehicle Camera", VisibleDefaultsOnly, BlueprintReadOnly, meta = (AllowPrivateAccess = "true"))
	class UCameraComponent* Camera;
};
