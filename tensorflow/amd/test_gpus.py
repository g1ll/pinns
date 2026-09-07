import os

# Configura o dispositivo AMD antes da inicialização do TensorFlow
os.environ.pop("CUDA_VISIBLE_DEVICES", None)
os.environ["HIP_VISIBLE_DEVICES"] = "0"
os.environ["HSA_OVERRIDE_GFX_VERSION"] = "10.3.0"
os.environ["LD_LIBRARY_PATH"] = "/opt/rocm/lib:" + os.environ.get("LD_LIBRARY_PATH", "")

import tensorflow as tf

# Validação do dispositivo
gpus = tf.config.list_physical_devices('GPU')
print("Dispositivos GPU detectados pelo TensorFlow:", gpus)
if gpus:
	tf.config.set_visible_devices(gpus[0], 'GPU')
	print("GPU AMD selecionada:", gpus[0].name)