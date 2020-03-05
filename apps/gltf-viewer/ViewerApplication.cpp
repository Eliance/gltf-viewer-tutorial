#include "ViewerApplication.hpp"

#include <iostream>
#include <numeric>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/io.hpp>
#include <glm/glm.hpp>

#include "utils/cameras.hpp"
#include "utils/gltf.hpp"
#include "utils/images.hpp"

#include <stb_image_write.h>
#include <tiny_gltf.h>

#define VERTEX_ATTRIB_POSITION_IDX  0
#define VERTEX_ATTRIB_NORMAL_IDX  1
#define VERTEX_ATTRIB_TEXCOORD0_IDX  2

#define ZERO 0.01f

void keyCallback(
    GLFWwindow *window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE) {
    glfwSetWindowShouldClose(window, 1);
  }
}

int ViewerApplication::run()
{
	//////////////////////////////////////////////////
	/// loading the scene
	//////////////////////////////////////////////////
	tinygltf::Model model;
	if (!loadGltfFile(model)){
		return -1;
	}

	//////////////////////////////////////////////////
	/// computing bounding box of the scene & camera vectors
	//////////////////////////////////////////////////
	glm::vec3 bboxMin;
	glm::vec3 bboxMax;
	computeSceneBounds(model, bboxMin, bboxMax);

	glm::vec3 center = glm::vec3((bboxMin.x+bboxMax.x)/2, (bboxMin.y+bboxMax.y)/2, (bboxMin.z+bboxMax.z)/2);
	glm::vec3 diag = glm::vec3(bboxMax.x-bboxMin.x, bboxMax.y-bboxMin.y, bboxMax.z-bboxMin.z);
	glm::vec3 eye; 
	//if a scene is flat:
	if ((bboxMax.z - bboxMin.z) < ZERO){
		eye = center + 2.f * glm::cross(diag, glm::vec3(0, 1, 0));
	} else {
		eye = center + diag;
	}

	float maxDistance = glm::length2(diag);
	maxDistance = maxDistance > ZERO ? maxDistance : 100.f;


	//////////////////////////////////////////////////
	/// camera Controller
	//////////////////////////////////////////////////
	
	std::unique_ptr<CameraController> cameraController = std::make_unique<TrackballCameraController>(m_GLFWHandle.window(), 0.5f * maxDistance);

	if (m_hasUserCamera) {
		cameraController->setCamera(m_userCamera);
	} else {
		const auto center = 0.5f * (bboxMax + bboxMin);
		const auto up = glm::vec3(0, 1, 0);
		const auto eye =
			diag.z > 0 ? center + diag : center + 2.f * glm::cross(diag, up);
		cameraController->setCamera(Camera{eye, center, up});
	}

	//////////////////////////////////////////////////
	// Loader shaders
	//////////////////////////////////////////////////
	const auto glslProgram =
		compileProgram({m_ShadersRootPath / m_AppName / m_vertexShader,
			m_ShadersRootPath / m_AppName / m_fragmentShader});

	//////////////////////////////////////////////////
	// compute textures
	//////////////////////////////////////////////////
	const std::vector<GLuint> textureObjects = createTextureObjects(model);
	GLuint whiteTexture = createWhiteTexture();

	//////////////////////////////////////////////////
	// compute matrix
	//////////////////////////////////////////////////

	const auto bufferObjects = createBufferObjects(model);
	std::vector<VaoRange> meshIndexToVaoRange;
	std::vector<GLuint> VAOs = createVertexArrayObjects(model, bufferObjects, meshIndexToVaoRange);

	//////////////////////////////////////////////////
	// uniform variables
	//////////////////////////////////////////////////

	// get uniform locations for uniform variables //
	
	const auto modelViewProjMatrixLocation =
		glGetUniformLocation(glslProgram.glId(), "uModelViewProjMatrix");
	const auto modelViewMatrixLocation =
		glGetUniformLocation(glslProgram.glId(), "uModelViewMatrix");
	const auto normalMatrixLocation =
		glGetUniformLocation(glslProgram.glId(), "uNormalMatrix");

	const auto uLightIntensityLocation =
		glGetUniformLocation(glslProgram.glId(), "uLightIntensity");
	const auto uLightDirectionLocation =
		glGetUniformLocation(glslProgram.glId(), "uLightDirection");

	const auto uBaseColorTextureLocation =
		glGetUniformLocation(glslProgram.glId(), "uBaseColorTexture");
		
	const auto uBaseColorFactorLocation =
		glGetUniformLocation(glslProgram.glId(), "uBaseColorFactor");
	if (uLightDirectionLocation < 0)
		std::cout << "no uLightDirection passed to fragment shader" << std::endl;
	if (uLightIntensityLocation < 0)
		std::cout << "no uLightIntensity passed to fragment shader" << std::endl;

	if (uBaseColorTextureLocation < 0)
		std::cout << "no uBaseColorTexture passed to fragment shader" << std::endl;

	const auto uMetallicFactorLocation = glGetUniformLocation(glslProgram.glId(), "uMetallicFactor");
	const auto uRoughnessFactorLocation = glGetUniformLocation(glslProgram.glId(), "uRoughnessFactor");
	const auto uMetallicRoughnessTextureLocation = glGetUniformLocation(glslProgram.glId(), "uMetallicRoughnessTexture");
	
	const auto uEmissiveTextureLocation = glGetUniformLocation(glslProgram.glId(), "uEmissiveTexture");
	const auto uEmissiveFactorLocation = glGetUniformLocation(glslProgram.glId(), "uEmissiveFactor");

	// declaration of matrix & vecs (light params)
	const auto projMatrix =
		glm::perspective(70.f, float(m_nWindowWidth) / m_nWindowHeight,
			0.001f * maxDistance, 1.5f * maxDistance);

	auto lightDirection = glm::vec3(0.f, 0.f, 1.f);
	auto lightIntensity = glm::vec3(1.f, 1.f, 1.f);
	bool lightFromCamera = false;
	
	// Setup OpenGL state for rendering
	glEnable(GL_DEPTH_TEST);
	glslProgram.use();

	//lambda function to bind textures
	

	
	const auto bindMaterial = [&](const auto materialIndex) {
		if (materialIndex >= 0){
			// only valid is materialIndex >= 0
			const auto &material = model.materials[materialIndex];
			const auto &pbrMetallicRoughness = material.pbrMetallicRoughness;

			if (uBaseColorTextureLocation >= 0) {
				auto textureObject = whiteTexture;
				const auto baseColorTextureIndex = pbrMetallicRoughness.baseColorTexture.index;
				if (baseColorTextureIndex >= 0) {
					// only valid if pbrMetallicRoughness.baseColorTexture.index >= 0:
					const tinygltf::Texture &texture = model.textures[baseColorTextureIndex];
					if (texture.source >= 0) {
						textureObject = textureObjects[texture.source];
					}
				}
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, textureObject);
				glUniform1i(uBaseColorTextureLocation, 0);
			}

			if (uBaseColorFactorLocation >= 0) {
				glUniform4f(uBaseColorFactorLocation,
					(float)pbrMetallicRoughness.baseColorFactor[0],
					(float)pbrMetallicRoughness.baseColorFactor[1],
					(float)pbrMetallicRoughness.baseColorFactor[2],
					(float)pbrMetallicRoughness.baseColorFactor[3]);
			}

			if (uMetallicFactorLocation >= 0) {
				glUniform1f(
					uMetallicFactorLocation, (float)pbrMetallicRoughness.metallicFactor);
			}
			if (uRoughnessFactorLocation >= 0) {
				glUniform1f(
					uRoughnessFactorLocation, (float)pbrMetallicRoughness.roughnessFactor);
			}
			if (uMetallicRoughnessTextureLocation > 0) {
				auto textureObject = 0u;
				if (pbrMetallicRoughness.metallicRoughnessTexture.index >= 0) {
					const auto &texture =
						model.textures[pbrMetallicRoughness.metallicRoughnessTexture
											.index];
					if (texture.source >= 0) {
						textureObject = textureObjects[texture.source];
					}
				}

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, textureObject);
				glUniform1i(uMetallicRoughnessTextureLocation, 1);
			}

			if (uEmissiveTextureLocation >= 0){

				auto textureObject = textureObjects[material.emissiveTexture.index];
				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, textureObject);
				glUniform1i(uEmissiveTextureLocation, 2);
			}
			if (uEmissiveFactorLocation >= 0){
				auto emissiveFactor = material.emissiveFactor;
				glUniform3f(uEmissiveFactorLocation,
					(float)emissiveFactor[0],
					(float)emissiveFactor[1],
					(float)emissiveFactor[2]
				);
			}
		} else {
			if (uBaseColorTextureLocation >= 0) {
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, whiteTexture);
				glUniform1i(uBaseColorTextureLocation, 0);
			}
			if (uBaseColorFactorLocation >= 0) {
				glUniform4f(uBaseColorFactorLocation, 1, 1, 1, 1);
			}
			if (uMetallicFactorLocation >= 0) {
				glUniform1f(uMetallicFactorLocation, 0.f);
			}
			if (uRoughnessFactorLocation >= 0) {
				glUniform1f(uRoughnessFactorLocation, 0.f);
			}
			if (uMetallicRoughnessTextureLocation >= 0) {
				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, whiteTexture);
				glUniform1i(uMetallicRoughnessTextureLocation, 0);
			}
			if (uEmissiveTextureLocation >= 0){
				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, whiteTexture);
				glUniform1i(uMetallicRoughnessTextureLocation, 0);
			}
			if (uMetallicFactorLocation >= 0){
				glUniform3f(uMetallicRoughnessTextureLocation, 0, 0, 0);
			}
		}
	};
	

	// Lambda function to draw the scene
  	const auto drawScene = [&](const Camera &camera) {
		glViewport(0, 0, m_nWindowWidth, m_nWindowHeight);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		const auto viewMatrix = camera.getViewMatrix();

		// The recursive function that should draw a node
		// We use a std::function because a simple lambda cannot be recursive
		const std::function<void(int, const glm::mat4 &)> drawNode = [&](int nodeIdx, const glm::mat4 &parentMatrix) {
			const auto &node = model.nodes[nodeIdx];
			glm::mat4 modelMatrix = getLocalToWorldMatrix(node, parentMatrix);
			int meshIdx = node.mesh;
			if (meshIdx >= 0){ //check if the node is a mesh
				glm::mat4 MV = viewMatrix * modelMatrix;
				glm::mat4 MVP = projMatrix*MV;
				glm::mat4 normalMatrix = transpose(inverse(MV));

				const auto lightDirectionNormalized = glm::normalize(
					glm::vec3(viewMatrix*glm::vec4(lightDirection, 0.f)));
				
				//sending uniform variables
				glUniformMatrix4fv(modelViewMatrixLocation, 1, GL_FALSE, glm::value_ptr(MV));
				glUniformMatrix4fv(modelViewProjMatrixLocation, 1, GL_FALSE, glm::value_ptr(MVP));
				glUniformMatrix4fv(normalMatrixLocation, 1, GL_FALSE, glm::value_ptr(normalMatrix));

				if (uLightIntensityLocation >= 0){
					glUniform3fv(uLightIntensityLocation, 1, glm::value_ptr(lightIntensity));
				}
				if (uLightDirectionLocation >= 0){
					if (lightFromCamera)
						glUniform3f(uLightDirectionLocation, 0, 0, 1);
					else
						glUniform3fv(uLightDirectionLocation, 1, glm::value_ptr(lightDirectionNormalized));
				}

				//material binding following primitives index.... todo
				if (uBaseColorTextureLocation >= 0){
					//glUniform1i();
				}

				//

				const auto &mesh = model.meshes[meshIdx];
				VaoRange vaoRange = meshIndexToVaoRange[meshIdx];
				const size_t primitiveNumber = mesh.primitives.size();
				for (size_t i=0; i < primitiveNumber; ++i){
					// Material binding
					bindMaterial(mesh.primitives[i].material);

					//vao binding
					const auto &vao = VAOs[vaoRange.begin + i];
					glBindVertexArray(vao);
					const auto &currentPrimitive = mesh.primitives[i];
					if (currentPrimitive.indices >= 0){
						const auto &accessor = model.accessors[currentPrimitive.indices];
						const auto &bufferView = model.bufferViews[accessor.bufferView];
						const auto byteOffset = accessor.byteOffset + bufferView.byteOffset;
						glDrawElements(currentPrimitive.mode, accessor.count, accessor.componentType, (const GLvoid *)byteOffset);
					} else {
						const auto accessorIdx = (*begin(currentPrimitive.attributes)).second;
						const auto &accessor = model.accessors[accessorIdx];
						glDrawArrays(currentPrimitive.mode, 0, accessor.count);
					}
				}
			}
			for (const auto &nodeIdx : node.children){
				drawNode(nodeIdx, modelMatrix);
			}
		};

		// Draw the scene referenced by gltf file
		if (model.defaultScene >= 0) {
			for (const auto &nodeIdx : model.scenes[model.defaultScene].nodes){
				drawNode(nodeIdx, glm::mat4(1));
			}
		}
	};


	if (!m_OutputPath.empty()){
		std::vector<unsigned char> pixels(m_nWindowWidth*m_nWindowHeight*3, 0);
		renderToImage(m_nWindowWidth, m_nWindowHeight, 3, pixels.data(), [&]() {drawScene(cameraController->getCamera());});
		flipImageYAxis(m_nWindowWidth, m_nWindowHeight, 3, pixels.data());
		const auto strPath = m_OutputPath.string();
		stbi_write_png(strPath.c_str(), m_nWindowWidth, m_nWindowHeight, 3, pixels.data(), 0);
		return 0;
	}

	// Loop until the user closes the window
	for (auto iterationCount = 0u; !m_GLFWHandle.shouldClose();
		++iterationCount) {
		const auto seconds = glfwGetTime();

		const auto camera = cameraController->getCamera();
		drawScene(camera);

		// GUI code:
		imguiNewFrame();
		{
			ImGui::Begin("GUI");
			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
				1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
			if (ImGui::CollapsingHeader("Camera", ImGuiTreeNodeFlags_DefaultOpen)) {
				ImGui::Text("eye: %.3f %.3f %.3f", camera.eye().x, camera.eye().y,
					camera.eye().z);
				ImGui::Text("center: %.3f %.3f %.3f", camera.center().x,
					camera.center().y, camera.center().z);
				ImGui::Text(
					"up: %.3f %.3f %.3f", camera.up().x, camera.up().y, camera.up().z);

				ImGui::Text("front: %.3f %.3f %.3f", camera.front().x, camera.front().y,
					camera.front().z);
				ImGui::Text("left: %.3f %.3f %.3f", camera.left().x, camera.left().y,
					camera.left().z);

				if (ImGui::Button("CLI camera args to clipboard")) {
				std::stringstream ss;
				ss << "--lookat " << camera.eye().x << "," << camera.eye().y << ","
					<< camera.eye().z << "," << camera.center().x << ","
					<< camera.center().y << "," << camera.center().z << ","
					<< camera.up().x << "," << camera.up().y << "," << camera.up().z;
				const auto str = ss.str();
				glfwSetClipboardString(m_GLFWHandle.window(), str.c_str());
				}

				static int cameraControllerType = 0;
				const auto cameraControllerTypeChanged =
					ImGui::RadioButton("Trackball", &cameraControllerType, 0) ||
					ImGui::RadioButton("First Person", &cameraControllerType, 1);
				if (cameraControllerTypeChanged) {
				const auto currentCamera = cameraController->getCamera();
				if (cameraControllerType == 0) {
					cameraController = std::make_unique<TrackballCameraController>(
						m_GLFWHandle.window(), 0.5f * maxDistance);
				} else {
					cameraController = std::make_unique<FirstPersonCameraController>(
						m_GLFWHandle.window(), 0.5f * maxDistance);
				}
				cameraController->setCamera(currentCamera);
				}
			}
			if (ImGui::CollapsingHeader("Light", ImGuiTreeNodeFlags_DefaultOpen)) {
				
				//light color & intensity
				static glm::vec3 lightColor(1.f, 1.f, 1.f);
				static float lightIntensityFactor = 1.f;
				if (ImGui::ColorEdit3("color", (float *)&lightColor) ||
					ImGui::SliderFloat("intensity", &lightIntensityFactor, 0.f, 10.f)) {
					lightIntensity = lightColor * lightIntensityFactor;
				}

				ImGui::Checkbox("light form camera", &lightFromCamera);

				if (!lightFromCamera){
					//light angles
					static float lightTheta = 0.f;
					static float lightPhi = 0.f;
					if (ImGui::SliderFloat("theta", &lightTheta, 0, glm::pi<float>()) ||
						ImGui::SliderFloat("phi", &lightPhi, 0, 2.f * glm::pi<float>())) {
						const auto sinPhi = glm::sin(lightPhi);
						const auto cosPhi = glm::cos(lightPhi);
						const auto sinTheta = glm::sin(lightTheta);
						const auto cosTheta = glm::cos(lightTheta);
						lightDirection = glm::vec3(sinTheta * cosPhi, cosTheta, sinTheta * sinPhi);
					}

				}
			}

			ImGui::End();
		}

		imguiRenderFrame();

		glfwPollEvents(); // Poll for and process events

		auto ellapsedTime = glfwGetTime() - seconds;
		auto guiHasFocus =
			ImGui::GetIO().WantCaptureMouse || ImGui::GetIO().WantCaptureKeyboard;
		if (!guiHasFocus) {
		cameraController->update(float(ellapsedTime));
		}

		m_GLFWHandle.swapBuffers(); // Swap front and back buffers
	}

	// TODO clean up allocated GL data

	return 0;
}

ViewerApplication::ViewerApplication(const fs::path &appPath, uint32_t width,
    uint32_t height, const fs::path &gltfFile,
    const std::vector<float> &lookatArgs, const std::string &vertexShader,
    const std::string &fragmentShader, const fs::path &output) :
    m_nWindowWidth(width),
    m_nWindowHeight(height),
    m_AppPath{appPath},
    m_AppName{m_AppPath.stem().string()},
    m_ImGuiIniFilename{m_AppName + ".imgui.ini"},
    m_ShadersRootPath{m_AppPath.parent_path() / "shaders"},
    m_gltfFilePath{gltfFile},
    m_OutputPath{output}
{
  if (!lookatArgs.empty()) {
    m_hasUserCamera = true;
    m_userCamera =
        Camera{glm::vec3(lookatArgs[0], lookatArgs[1], lookatArgs[2]),
            glm::vec3(lookatArgs[3], lookatArgs[4], lookatArgs[5]),
            glm::vec3(lookatArgs[6], lookatArgs[7], lookatArgs[8])};
  }

  if (!vertexShader.empty()) {
    m_vertexShader = vertexShader;
  }

  if (!fragmentShader.empty()) {
    m_fragmentShader = fragmentShader;
  }

  ImGui::GetIO().IniFilename =
      m_ImGuiIniFilename.c_str(); // At exit, ImGUI will store its windows
                                  // positions in this file

  glfwSetKeyCallback(m_GLFWHandle.window(), keyCallback);

  printGLVersion();
}

bool ViewerApplication::loadGltfFile(tinygltf::Model & model){
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;
    bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, m_gltfFilePath.string());

    if (!warn.empty()) {
        std::cout << "Warn: " << warn << std::endl;
    }
    if (!err.empty()) {
        std::cerr << err << std::endl;
    }

    if (!ret) {
		std::cerr << "Failed to parse glTF" << std::endl;
		return false;
    }

    return ret;
}

std::vector<GLuint> ViewerApplication::createBufferObjects(const tinygltf::Model &model){
    size_t size = model.buffers.size();
    std::vector<GLuint> bufferObjects(size, 0);
    glGenBuffers(size, bufferObjects.data());
    
    for (size_t i=0; i < size; ++i){
        glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[i]);
        glBufferStorage(GL_ARRAY_BUFFER, model.buffers[i].data.size(), model.buffers[i].data.data(), 0);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0); 

    return bufferObjects;
}

std::vector<GLuint> ViewerApplication::createVertexArrayObjects( const tinygltf::Model &model, const std::vector<GLuint> &bufferObjects, std::vector<ViewerApplication::VaoRange> & meshIndexToVaoRange){
    std::vector<GLuint> vertexArrayObjects;

    for (const auto &mesh : model.meshes){
		GLsizei primitiveNumber = mesh.primitives.size();
		GLsizei offset = vertexArrayObjects.size();
		vertexArrayObjects.resize(offset + primitiveNumber);
		meshIndexToVaoRange.push_back(ViewerApplication::VaoRange{offset, primitiveNumber});
		
		glGenVertexArrays(primitiveNumber ,vertexArrayObjects.data()+offset);

		for (size_t i=0; i < primitiveNumber; ++i){
			glBindVertexArray(vertexArrayObjects[offset + i]);
			auto currentPrimitive = mesh.primitives[i];

			{ // I'm opening a scope because I want to reuse the variable iterator in the code for NORMAL and TEXCOORD_0
				const auto iterator = currentPrimitive.attributes.find("POSITION");
				if (iterator != end(currentPrimitive.attributes)) { // If "POSITION" has been found in the map
					const auto accessorIdx = (*iterator).second; // (*iterator).first is the key "POSITION", (*iterator).second is the value, ie. the index of the accessor for this attribute
					const auto &accessor = model.accessors[accessorIdx]; //get the correct tinygltf::Accessor from model.accessors
					const auto &bufferView = model.bufferViews[accessor.bufferView]; //get the correct tinygltf::BufferView from model.bufferViews. You need to use the accessor
					const int bufferIdx = bufferView.buffer; //get the index of the buffer used by the bufferView (you need to use it)
					//const int bufferIdx = model.buffers[bufferView.buffer]; //get the index of the buffer used by the bufferView (you need to use it)
					
					const auto bufferObject = bufferObjects[bufferIdx]; //get the correct buffer object from the buffer index

					glEnableVertexAttribArray(VERTEX_ATTRIB_POSITION_IDX); //Enable the vertex attrib array corresponding to POSITION with glEnableVertexAttribArray (you need to use VERTEX_ATTRIB_POSITION_IDX which is defined at the top of the file)
					glBindBuffer(GL_ARRAY_BUFFER, bufferObject); //Bind the buffer object to GL_ARRAY_BUFFER

					const auto byteOffset = bufferView.byteOffset + accessor.byteOffset; //Compute the total byte offset using the accessor and the buffer view
					glVertexAttribPointer(VERTEX_ATTRIB_POSITION_IDX, accessor.type, GL_FLOAT, GL_FALSE, GLsizei(bufferView.byteStride), (const GLvoid *)byteOffset);// TODO Call glVertexAttribPointer with the correct arguments. 
					// Remember size is obtained with accessor.type, type is obtained with accessor.componentType. 
					// The stride is obtained in the bufferView, normalized is always GL_FALSE, and pointer is the byteOffset (don't forget the cast).
				}
			}
			{
				const auto iterator = currentPrimitive.attributes.find("NORMAL");
				if (iterator != end(currentPrimitive.attributes)) { 
					const auto accessorIdx = (*iterator).second; // (*iterator).first is the key "POSITION", (*iterator).second is the value, ie. the index of the accessor for this attribute
					const auto &accessor = model.accessors[accessorIdx]; //get the correct tinygltf::Accessor from model.accessors
					const auto &bufferView = model.bufferViews[accessor.bufferView]; //get the correct tinygltf::BufferView from model.bufferViews. You need to use the accessor
					const int bufferIdx = bufferView.buffer; //get the index of the buffer used by the bufferView (you need to use it)
					//const int bufferIdx = model.buffers[bufferView.buffer]; //get the index of the buffer used by the bufferView (you need to use it)
					
					const auto bufferObject = bufferObjects[bufferIdx]; //get the correct buffer object from the buffer index

					glEnableVertexAttribArray(VERTEX_ATTRIB_NORMAL_IDX); //Enable the vertex attrib array corresponding to POSITION with glEnableVertexAttribArray (you need to use VERTEX_ATTRIB_POSITION_IDX which is defined at the top of the file)
					glBindBuffer(GL_ARRAY_BUFFER, bufferObject); //Bind the buffer object to GL_ARRAY_BUFFER

					const auto byteOffset = bufferView.byteOffset + accessor.byteOffset; //Compute the total byte offset using the accessor and the buffer view
					glVertexAttribPointer(VERTEX_ATTRIB_NORMAL_IDX, accessor.type, GL_FLOAT, GL_FALSE, GLsizei(bufferView.byteStride), (const GLvoid *)byteOffset);// TODO Call glVertexAttribPointer with the correct arguments. 
					// Remember size is obtained with accessor.type, type is obtained with accessor.componentType. 
					// The stride is obtained in the bufferView, normalized is always GL_FALSE, and pointer is the byteOffset (don't forget the cast).
				}
			}
			{ // TEXCOORD_0 attribute
					const auto iterator = currentPrimitive.attributes.find("TEXCOORD_0");
					if (iterator != end(currentPrimitive.attributes)) {
					const auto accessorIdx = (*iterator).second;
					const auto &accessor = model.accessors[accessorIdx];
					const auto &bufferView = model.bufferViews[accessor.bufferView];
					const auto bufferIdx = bufferView.buffer;

					glEnableVertexAttribArray(VERTEX_ATTRIB_TEXCOORD0_IDX);
					assert(GL_ARRAY_BUFFER == bufferView.target);
					glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[bufferIdx]);
					glVertexAttribPointer(VERTEX_ATTRIB_TEXCOORD0_IDX, accessor.type,
						accessor.componentType, GL_FALSE, GLsizei(bufferView.byteStride),
						(const GLvoid *)(accessor.byteOffset + bufferView.byteOffset));
			}
			}
			// Index array if defined (if the primitive is use several times)
			uint primitiveIndices = currentPrimitive.indices;
			if (primitiveIndices >= 0) {
				const auto accessorIdx = primitiveIndices;
				const auto &accessor = model.accessors[accessorIdx];
				const auto &bufferView = model.bufferViews[accessor.bufferView];
				const auto bufferIdx = bufferView.buffer;

				const auto bufferObject = bufferObjects[bufferIdx];

				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferObject);
			}
		}

    }
	glBindVertexArray(0);
	//std::clog << "Number of VAOs: " << vertexArrayObjects.size() << std::endl;
    return vertexArrayObjects;
}



std::vector<GLuint> ViewerApplication::createTextureObjects(const tinygltf::Model &model) const{

	const auto & modelTextures = model.textures;
	size_t nbTexturesToGen = modelTextures.size();
	
	// Generate the texture object
	auto texObject = std::vector<GLuint>(nbTexturesToGen, 0);

	tinygltf::Sampler defaultSampler;
	defaultSampler.minFilter = GL_LINEAR;
	defaultSampler.magFilter = GL_LINEAR;
	defaultSampler.wrapS = GL_REPEAT;
	defaultSampler.wrapT = GL_REPEAT;
	defaultSampler.wrapR = GL_REPEAT;

	glActiveTexture(GL_TEXTURE0);

	glGenTextures(nbTexturesToGen, texObject.data());

	//Texture from the 3D model passed in argument :
	for (size_t i = 0; i < nbTexturesToGen; ++i){
		glBindTexture(GL_TEXTURE_2D, texObject[i]);
			
		/*
		// Set sampling parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		// Set wrapping parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
		*/

		const auto &texture = modelTextures[i]; // get i-th texture
		assert(texture.source >= 0); // ensure a source image is present
		const auto &image = model.images[texture.source]; // get the image
		const auto &sampler = texture.sampler >= 0 ? model.samplers[texture.sampler] : defaultSampler;
		// fill the texture object with the data from the image
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width, image.height, 0, GL_RGBA, image.pixel_type, image.image.data());
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
		sampler.minFilter != -1 ? sampler.minFilter : GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
			sampler.magFilter != -1 ? sampler.magFilter : GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, sampler.wrapS);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, sampler.wrapT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, sampler.wrapR);

		if (sampler.minFilter == GL_NEAREST_MIPMAP_NEAREST ||
			sampler.minFilter == GL_NEAREST_MIPMAP_LINEAR ||
			sampler.minFilter == GL_LINEAR_MIPMAP_NEAREST ||
			sampler.minFilter == GL_LINEAR_MIPMAP_LINEAR) {
			glGenerateMipmap(GL_TEXTURE_2D);
		}
	}
	glBindTexture(GL_TEXTURE_2D, 0);

	return texObject;
}

GLuint ViewerApplication::createWhiteTexture() const{
	GLuint texObject;

	// Set sampling parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	// Set wrapping parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
	
	glGenTextures(1, &texObject);
	glBindTexture(GL_TEXTURE_2D, texObject); // Bind to target GL_TEXTURE_2D
		float white[] = {1, 1, 1, 1}; // "Useless" black image of RGB (3 components) float values
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 2, 2, 0, GL_RGB, GL_FLOAT, white); // Set image data
	glBindTexture(GL_TEXTURE_2D, 0);

	return texObject;
}