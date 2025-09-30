# CFD_SMOKE2D_DIVMLS
⚠️ 본 저장소에는 지도 교수님께서 제공하신 유체 시뮬레이션 관련 코드가 비공개로 유지되어, 연구에서 사용된 전체 코드가 포함되어 있지 않습니다.
이로 인해 본 저장소의 코드만으로는 프로젝트를 실행할 수 없고, DIVMLS와 GPU 연산 및 관련 서브모듈만 포함되어 있음을 알려드립니다.
<hr>

# 🗃️ Summary 
Semi-Lagrangian 이류 과정에서 역추적(Backward tracing)한 위치의 주변 속도를 **Divergence-constrained MLS(Moving least squares)를 이용하여 보간하고 그 결과를 이류된 속도 데이터의 외력으로 적용**해 연기 시뮬레이션의 난류 표현을 개선한다.

<img width="578" height="667" alt="image" src="https://github.com/user-attachments/assets/56d75a11-97d9-4aad-9124-0651853c8e3b" />

Figure. Simulate smoke with two colliding densities. (a) Stable Fluid, (b)DCMLS interpolation, (c) DCMLS-applied external force

- Stable Fluids의 안정적인 경로를 유지하며 난류 강화
- 고차보간으로 인해 발생하는 Noise 완화
- 특이값 분해(SVD)를 병렬 프로그래밍에 적합한 알고리즘으로 대체한다면 충분히 고속화 가능
<br>

# 연구 목표
## 💡연구의 필요성
<img width="623" height="346" alt="image" src="https://github.com/user-attachments/assets/3a9b4f9f-2001-423c-bd57-1a940f614835" />

Figure. Comparison of turbulence generated using the previous MLS interpolation : (a) Stable Fluids[1], (b) DCMLS interpolation[3].

**[Stable Fluids 기법의 단점]** <br>
유체 시뮬레이션에 주로 사용되는 Stable Fluids 기법은 Fluids기법은 세미-라그랑주(Semi-Lagrangian) 이류를 사용해 안정적인 시뮬레이션 결과를 보장함
그러나, **수치 손실로 인해 디테일한 난류 표현이 어려움**

**[DCMLS로 보간법을 대체한 기법의 단점]** <br>
이를 개선하기 위해 DCMLS(Divergence-constrained MLS)로 기존 이류 단계의 보간법을 대체해 난류를 강화한 연구가 있지만, <br>
결과가 **매우 노이즈하고 불안정**할 뿐만 아니라 일부 영역에서 튀는 현상이 발생하는 단점이 있음

# 제안 기법

### DCMLS란?

**Divergence-free(발산이 0인 상태, 물질의 양이 변하지 않는 것)를** 고려해 벡터 충돌 시, 서로 상쇄되지 않고 운동량이 보존되도록 한 보간법
<br>
`DCMLS로 보간을 대체한 난류가 노이즈하고 불안정한 원인` : 난류가 유체의 전체적으로 강화되기 때문

## ✅ 적용 방법

<img width="1297" height="251" alt="image" src="https://github.com/user-attachments/assets/936e08e8-87aa-4009-9972-234d9aa2bb2a" />

1. 이중선형보간법과 DCMLS로 보간된 속도장을 각각 연산
2. 이중선형보간법 기반 속도장에 DCMLS 기반 속도장을 외력으로 적용

# 결과
