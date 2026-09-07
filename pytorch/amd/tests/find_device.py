import torch


def main():
    print(f"PyTorch: {torch.__version__}")
    print(f"HIP/ROCm: {torch.version.hip}")

    if torch.version.hip is None:
        raise RuntimeError(
            "A instalacao atual nao e uma build ROCm/HIP; CUDA NVIDIA nao e aceita."
        )

    if not torch.cuda.is_available():
        raise RuntimeError("Nenhuma GPU AMD foi detectada pelo PyTorch/ROCm.")

    if torch.cuda.device_count() < 1:
        raise RuntimeError("O PyTorch/ROCm nao encontrou nenhum dispositivo GPU.")

    device = torch.device("cuda:0")
    device_name = torch.cuda.get_device_name(device)
    device_name_lower = device_name.lower()
    if "amd" not in device_name_lower and "radeon" not in device_name_lower:
        raise RuntimeError(
            f"Dispositivo inesperado: {device_name}. Esperado: GPU AMD Radeon."
        )

    values = torch.tensor([1.0, 2.0, 3.0], device=device)
    result = values.square()
    torch.cuda.synchronize(device)

    if result.device != device:
        raise RuntimeError(
            f"A operacao nao foi executada na GPU: resultado em {result.device}."
        )

    print(f"GPU AMD: {device_name}")
    print(f"Dispositivo usado: {result.device}")
    print(f"Operacao na GPU: {result.tolist()}")
    print("Teste ROCm AMD concluido com sucesso.")


if __name__ == "__main__":
    main()