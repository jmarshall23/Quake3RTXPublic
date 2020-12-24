// d3d12_raytrace_mega.cpp
//

#include "d3d12_local.h"
#include "ImagePacker.h"

struct MegaTexture_t {
	ComPtr<ID3D12Resource> textureUploadHeap;
	ComPtr<ID3D12Resource> texture2D;
};

MegaTexture_t megaTexture;
MegaTexture_t megaTextureNormal;

const int r_megaTextureSize = 16384;

extern ComPtr<ID3D12DescriptorHeap> m_srvUavHeap;

iceMegaTexture* diffuseMegaTexture = NULL;
iceMegaTexture* normalMegaTexture = NULL;

void GL_InitMegaTextures(void) {
	if (diffuseMegaTexture != NULL) {
		delete diffuseMegaTexture;
		diffuseMegaTexture = NULL;
	}

	if (normalMegaTexture != NULL) {
		delete normalMegaTexture;
		normalMegaTexture = NULL;
	}

	diffuseMegaTexture = new iceMegaTexture();
	diffuseMegaTexture->InitTexture();

	normalMegaTexture = new iceMegaTexture();
	normalMegaTexture->InitTexture();
}

void GL_RegisterTexture(const char* texturePath, int width, int height, byte* data) {
	if (diffuseMegaTexture == NULL)
		return;

	diffuseMegaTexture->RegisterTexture(texturePath, width, height, data);
}

void GL_LoadMegaTexture(D3D12_CPU_DESCRIPTOR_HANDLE& srvPtr) {
	int width = r_megaTextureSize;
	int height = r_megaTextureSize;

	//byte* buffer;
	//LoadTGA("mega/mega.tga", &buffer, &width, &height, NULL);
	//FILE* f = fopen("base/mega/mega.raw", "rb");
	//fread(buffer, 1, width * height * 4, f);
	//fclose(f);

	diffuseMegaTexture->BuildMegaTexture();
	byte* buffer = diffuseMegaTexture->GetMegaBuffer();

	// Create the texture.
	{
		// Describe and create a Texture2D.
		D3D12_RESOURCE_DESC textureDesc = {};
		textureDesc.MipLevels = 1;
		textureDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
		textureDesc.Width = width;
		textureDesc.Height = height;
		textureDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
		textureDesc.DepthOrArraySize = 1;
		textureDesc.SampleDesc.Count = 1;
		textureDesc.SampleDesc.Quality = 0;
		textureDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;

		ThrowIfFailed(m_device->CreateCommittedResource(
			&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT),
			D3D12_HEAP_FLAG_NONE,
			&textureDesc,
			D3D12_RESOURCE_STATE_COPY_DEST,
			nullptr,
			IID_PPV_ARGS(&megaTexture.texture2D)));

		const UINT64 uploadBufferSize = GetRequiredIntermediateSize(megaTexture.texture2D.Get(), 0, 1);

		// Create the GPU upload buffer.
		ThrowIfFailed(m_device->CreateCommittedResource(
			&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD),
			D3D12_HEAP_FLAG_NONE,
			&CD3DX12_RESOURCE_DESC::Buffer(uploadBufferSize),
			D3D12_RESOURCE_STATE_GENERIC_READ,
			nullptr,
			IID_PPV_ARGS(&megaTexture.textureUploadHeap)));

		D3D12_SUBRESOURCE_DATA textureData = {};
		textureData.pData = &buffer[0];
		textureData.RowPitch = width * 4;
		textureData.SlicePitch = textureData.RowPitch * height;

		UpdateSubresources(m_commandList.Get(), megaTexture.texture2D.Get(), megaTexture.textureUploadHeap.Get(), 0, 0, 1, &textureData);
		m_commandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(megaTexture.texture2D.Get(), D3D12_RESOURCE_STATE_COPY_DEST, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE));

		// Describe and create a SRV for the texture.
		D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
		srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
		srvDesc.Format = textureDesc.Format;
		srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
		srvDesc.Texture2D.MipLevels = 1;
		m_device->CreateShaderResourceView(megaTexture.texture2D.Get(), &srvDesc, srvPtr);
	}
}


void GL_LoadMegaNormalTexture(D3D12_CPU_DESCRIPTOR_HANDLE& srvPtr) {
	int width = r_megaTextureSize;
	int height = r_megaTextureSize;

	//byte* buffer;
	//LoadTGA("mega/mega_local.tga", &buffer, &width, &height, NULL);
	//FILE* f = fopen("base/mega/mega.raw", "rb");
	//fread(buffer, 1, width * height * 4, f);
	//fclose(f);

	normalMegaTexture->BuildMegaTexture();
	byte* buffer = normalMegaTexture->GetMegaBuffer();

	// Create the texture.
	{
		// Describe and create a Texture2D.
		D3D12_RESOURCE_DESC textureDesc = {};
		textureDesc.MipLevels = 1;
		textureDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
		textureDesc.Width = width;
		textureDesc.Height = height;
		textureDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
		textureDesc.DepthOrArraySize = 1;
		textureDesc.SampleDesc.Count = 1;
		textureDesc.SampleDesc.Quality = 0;
		textureDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;

		ThrowIfFailed(m_device->CreateCommittedResource(
			&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT),
			D3D12_HEAP_FLAG_NONE,
			&textureDesc,
			D3D12_RESOURCE_STATE_COPY_DEST,
			nullptr,
			IID_PPV_ARGS(&megaTextureNormal.texture2D)));

		const UINT64 uploadBufferSize = GetRequiredIntermediateSize(megaTextureNormal.texture2D.Get(), 0, 1);

		// Create the GPU upload buffer.
		ThrowIfFailed(m_device->CreateCommittedResource(
			&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD),
			D3D12_HEAP_FLAG_NONE,
			&CD3DX12_RESOURCE_DESC::Buffer(uploadBufferSize),
			D3D12_RESOURCE_STATE_GENERIC_READ,
			nullptr,
			IID_PPV_ARGS(&megaTextureNormal.textureUploadHeap)));

		D3D12_SUBRESOURCE_DATA textureData = {};
		textureData.pData = &buffer[0];
		textureData.RowPitch = width * 4;
		textureData.SlicePitch = textureData.RowPitch * height;

		UpdateSubresources(m_commandList.Get(), megaTextureNormal.texture2D.Get(), megaTextureNormal.textureUploadHeap.Get(), 0, 0, 1, &textureData);
		m_commandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(megaTextureNormal.texture2D.Get(), D3D12_RESOURCE_STATE_COPY_DEST, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE));

		// Describe and create a SRV for the texture.
		D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
		srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
		srvDesc.Format = textureDesc.Format;
		srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
		srvDesc.Texture2D.MipLevels = 1;
		m_device->CreateShaderResourceView(megaTextureNormal.texture2D.Get(), &srvDesc, srvPtr);
	}
}
/*
==============
GL_FindMegaTile
==============
*/
void GL_FindMegaTile(const char* name, float& x, float& y, float& width, float& height) {
	if (diffuseMegaTexture == NULL)
		return;

	diffuseMegaTexture->FindMegaTile(name, x, y, width, height);
}
/*
==============
GL_FindMegaTile
==============
*/
void GL_FindMegaTile(const char* name, float* x, float* y, float* width, float* height) {
	float _x, _y, _width, _height;
	if (diffuseMegaTexture == NULL)
		return;

	diffuseMegaTexture->FindMegaTile(name, _x, _y, _width, _height);

	*x = _x;
	*y = _y;
	*width = _width;
	*height = _height;
}

/*
==============
R_CopyImage
==============
*/
void R_CopyImage(byte* source, int sourceX, int sourceY, int sourceWidth, byte* dest, int destX, int destY, int destWidth, int width, int height)
{
	for (int y = 0; y < height; y++)
	{
		int _x = 0;
		int _y = y * 4;
		int destPos = (destWidth * (_y + (destY * 4))) + (_x + (destX * 4));
		int sourcePos = (sourceWidth * (_y + (sourceY * 4))) + (_x + (sourceX * 4));

		memcpy(&dest[destPos], &source[sourcePos], width * 4);
	}
}

/*
=====================
iceMegaTexture::iceMegaTexture
=====================
*/
iceMegaTexture::iceMegaTexture() {
	isRegistered = false;
	megaTextureBuffer = NULL;
	imagePacker = NULL;
}
/*
=====================
iceMegaTexture::iceMegaTexture
=====================
*/
iceMegaTexture::~iceMegaTexture() {
	for (int i = 0; i < megaEntries.size(); i++)
	{
		delete megaEntries[i].data_copy;
		megaEntries[i].data_copy = NULL;
	}
	megaEntries.clear();

	if (megaTextureBuffer != NULL) {
		delete megaTextureBuffer;
		megaTextureBuffer = NULL;
	}

	if (imagePacker != NULL) {
		delete imagePacker;
	}
}

/*
=======================
iceMegaTexture::InitTexture
=======================
*/
void iceMegaTexture::InitTexture(void) {
	int megaSize = r_megaTextureSize;

	Com_Printf("Init MegaTexture %dx%d\n", megaSize, megaSize);
	imagePacker = new idImagePacker(megaSize, megaSize);

	megaTextureBuffer = new byte[megaSize * megaSize * 4];
	memset(megaTextureBuffer, 0, megaSize * megaSize * 4);
}

/*
=======================
iceMegaTexture::RegisterTexture
=======================
*/
void iceMegaTexture::RegisterTexture(const char* texturePath, int width, int height, byte* data) {
	int tileId = megaEntries.size();

	if (isRegistered) {
		Com_Printf("iceMegaTexture::RegisterTexture: %s trying to be registered outside of registration!\n", texturePath);
		return;
	}

	//common->Printf("RegisterTexture: %s\n", texturePath);

	idSubImage subImage = imagePacker->PackImage(width, height, false);

	// Check to make sure we don't have any megaEntries
	iceMegaEntry newEntry;
	strcpy(newEntry.texturePath, texturePath);
	newEntry.width = width;
	newEntry.height = height;
	newEntry.x = subImage.x;
	newEntry.y = subImage.y;
	newEntry.data_copy = new byte[width * height * 4];
	memcpy(newEntry.data_copy, data, width * height * 4);

	megaEntries.push_back(newEntry);
}

/*
=======================
iceMegaTexture::BuildMegaTexture
=======================
*/
void iceMegaTexture::BuildMegaTexture(void) {
	int megaSize = r_megaTextureSize;

	// Update all of our megatexture entries.
	Com_Printf("Updating %d entries...\n", megaEntries.size());
	for (int i = 0; i < megaEntries.size(); i++)
	{
		if (megaEntries[i].x + megaEntries[i].width > r_megaTextureSize) {
			continue;
		}

		if (megaEntries[i].y + megaEntries[i].height > r_megaTextureSize) {
			continue;
		}
		R_CopyImage(megaEntries[i].data_copy, 0, 0, megaEntries[i].width, megaTextureBuffer, megaEntries[i].x, megaEntries[i].y, r_megaTextureSize, megaEntries[i].width, megaEntries[i].height);
	}

	// R_WriteTGA("testme.tga", megaTextureBuffer, r_megaTextureSize.GetInteger(), r_megaTextureSize.GetInteger());

	isRegistered = true;
}

/*
=======================
iceMegaTexture::FindMegaTile
=======================
*/
void iceMegaTexture::FindMegaTile(const char* name, float& x, float& y, float& width, float& height)
{
	for (int i = 0; i < megaEntries.size(); i++) {
		if (megaEntries[i].texturePath == name) {
			x = megaEntries[i].x;
			y = megaEntries[i].y;
			width = megaEntries[i].width;
			height = megaEntries[i].height;
			return;
		}
	}

	for (int i = 0; i < megaEntries.size(); i++) {
		if (strstr(megaEntries[i].texturePath, name)) {
			x = megaEntries[i].x;
			y = megaEntries[i].y;
			width = megaEntries[i].width;
			height = megaEntries[i].height;
			return;
		}
	}
	x = -1;
	y = -1;
	width = -1;
	height = -1;
}