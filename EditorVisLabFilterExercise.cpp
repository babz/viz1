// Copyright (C) 2011 vis-group. All Rights Reserved.


#include "EditorVislabFilterExercise.h"
#include "vrc_plugin_editor_vislabfilterexercise.cpp"
#include "svn_plugin_editor_vislabfilterexercise.cpp"
#include "CoreVersion.h"
#include "VolumeShop.h"
#include "Exception.h"
#include "Handle.h"
#include "moc_EditorVisLabFilterExercise.cpp"

DEFINE_EDITORPLUGIN(EditorVisLabFilterExercise)
DEFINE_PLUGIN_VERSION(EditorVisLabFilterExercise,0,0,1)

/**
 * gradient by central difference method
 **/
template <class dstType, unsigned int dstN>
class BackwardDifferenceGradientFilter
{
public:
	template <class srcType, unsigned int srcN>
	void transform(Voxel<dstType,dstN> & OutputVoxel, const Voxel<srcType,srcN> Neighborhood[3][3][3])
	{
		const Voxel<float,1> CurrentVoxel( Neighborhood[1][1][1].Get(3) );

		const Voxel<float,1> MinusVoxelX( Neighborhood[1][1][0].Get(3) );
		const Voxel<float,1> MinusVoxelY( Neighborhood[1][0][1].Get(3) );
		const Voxel<float,1> MinusVoxelZ( Neighborhood[0][1][1].Get(3) );

		// Backward Differences are used for gradient computation
		// http://en.wikipedia.org/wiki/Finite_difference
		OutputVoxel.Set(0,(CurrentVoxel-MinusVoxelX).Get(0) );
		OutputVoxel.Set(1,(CurrentVoxel-MinusVoxelY).Get(0) );
		OutputVoxel.Set(2,(CurrentVoxel-MinusVoxelZ).Get(0) );
	};
};

/**
 * gradient by central difference method
 **/
template <class dstType, unsigned int dstN>
class CentralDifferenceGradientFilter
{
public:
	template <class srcType, unsigned int srcN>
	void transform(Voxel<dstType,dstN> & OutputVoxel, const Voxel<srcType,srcN> Neighborhood[3][3][3])
	{
		const Voxel<float,1> CurrentVoxel( Neighborhood[1][1][1].Get(3) );

		const Voxel<float,1> MinusVoxelX( Neighborhood[1][1][0].Get(3) );
		const Voxel<float,1> PlusVoxelX( Neighborhood[1][1][2].Get(3) );

		const Voxel<float,1> MinusVoxelY( Neighborhood[1][0][1].Get(3) );
		const Voxel<float,1> PlusVoxelY( Neighborhood[1][2][1].Get(3) );

		const Voxel<float,1> MinusVoxelZ( Neighborhood[0][1][1].Get(3) );
		const Voxel<float,1> PlusVoxelZ( Neighborhood[2][1][1].Get(3) );

		// TODO 1: Implement Gradient by Central Differences
		// http://www.cs.utah.edu/~gk/MS/html/node33.html
		// currently Backward Differences are used for gradient computation
		// http://en.wikipedia.org/wiki/Finite_difference
		
		OutputVoxel.Set(0, 0.5f * (PlusVoxelX-MinusVoxelX).Get(0));
		OutputVoxel.Set(1, 0.5f * (PlusVoxelY-MinusVoxelY).Get(0));
		OutputVoxel.Set(2, 0.5f * (PlusVoxelZ-MinusVoxelZ).Get(0));

		// END TODO 1
	};
};

/**
 * gradient by sobel
 **/
template <class dstType, unsigned int dstN>
class SobelGradientFilter
{
public:
	template <class srcType, unsigned int srcN>
	void transform(Voxel<dstType,dstN> & OutputVoxel, const Voxel<srcType,srcN> Neighborhood[3][3][3])
	{
		// upper slice
		const Voxel<float,1> v220( Neighborhood[2][2][0].Get(3) ); const Voxel<float,1> v221( Neighborhood[2][2][1].Get(3) ); const Voxel<float,1> v222( Neighborhood[2][2][2].Get(3) );
		const Voxel<float,1> v210( Neighborhood[2][1][0].Get(3) ); const Voxel<float,1> v211( Neighborhood[2][1][1].Get(3) ); const Voxel<float,1> v212( Neighborhood[2][1][2].Get(3) );
		const Voxel<float,1> v200( Neighborhood[2][0][0].Get(3) ); const Voxel<float,1> v201( Neighborhood[2][0][1].Get(3) ); const Voxel<float,1> v202( Neighborhood[2][0][2].Get(3) );

		// current slice
		const Voxel<float,1> v120( Neighborhood[1][2][0].Get(3) ); const Voxel<float,1> v121( Neighborhood[1][2][1].Get(3) ); const Voxel<float,1> v122( Neighborhood[1][2][2].Get(3) );
		const Voxel<float,1> v110( Neighborhood[1][1][0].Get(3) ); const Voxel<float,1> v111( Neighborhood[1][1][1].Get(3) ); const Voxel<float,1> v112( Neighborhood[1][1][2].Get(3) );
		const Voxel<float,1> v100( Neighborhood[1][0][0].Get(3) ); const Voxel<float,1> v101( Neighborhood[1][0][1].Get(3) ); const Voxel<float,1> v102( Neighborhood[1][0][2].Get(3) );

		// lower slice
		const Voxel<float,1> v020( Neighborhood[0][2][0].Get(3) ); const Voxel<float,1> v021( Neighborhood[0][2][1].Get(3) ); const Voxel<float,1> v022( Neighborhood[0][2][2].Get(3) );
		const Voxel<float,1> v010( Neighborhood[0][1][0].Get(3) ); const Voxel<float,1> v011( Neighborhood[0][1][1].Get(3) ); const Voxel<float,1> v012( Neighborhood[0][1][2].Get(3) );
		const Voxel<float,1> v000( Neighborhood[0][0][0].Get(3) ); const Voxel<float,1> v001( Neighborhood[0][0][1].Get(3) ); const Voxel<float,1> v002( Neighborhood[0][0][2].Get(3) );

		// TODO 2: Gradient by Sobel Filter
		// Implement gradient computation using the Sobel operator
		// http://en.wikipedia.org/wiki/Sobel_operator
//z-Filterkernel: (x columns / y rows)
		// -1 -2 -1		0 0 0	1 2 1
		// -2 -4 -2		0 0 0	2 4 2
		// -1 -2 -1		0 0 0	1 2 1

		// x, y, z
		float filterZ[3][3][3] = { -1, -2, -1, -2, -4, -2, -1, -2, -1,
									0, 0, 0, 0, 0, 0, 0, 0, 0,
									1, 2, 1, 2, 4, 2, 1, 2, 1 };

		//y-Filterkernel (x columns / z rows)
		// -1 -2 -1		-2 -4 -2	-1 -2 -1
		//  0  0  0		 0  0  0	 0  0  0
		//  1  2  1		 2  4  2	 1  2  1

		float filterY[3][3][3] = { -1, -2, -1, 0, 0, 0, 1, 2, 1,
								   -2, -4, -2, 0, 0, 0, 2, 4, 2,
								   -1, -2, -1, 0, 0, 0, 1, 2, 1 };

		//x-Filterkernel (z columns / y rows)
		// -1 0 1	-2 0 2	 -1 0 1
		// -2 0 2	-4 0 4	 -2 0 2
		// -1 0 1	-2 0 2	 -1 0 1

		float filterX[3][3][3] = { -1, 0, 1, -2, 0, 2, -1, 0, 1,
								   -2, 0, 2, -4, 0, 4, -2, 0, 2,
								   -1, 0, 1, -2, 0, 2, -1, 0, 1 };

		float gZ = 0.0f; 
		float gY = 0.0f; 
		float gX = 0.0f;

		for (int x = 0; x < 3; ++x) {
			for (int y = 0; y < 3; ++y) {
				for (int z = 0; z < 3; ++z) {
					gZ += filterZ[z][y][x] * Neighborhood[z][y][x].Get(3);
				}
			}
		}
		for (int x = 0; x < 3; ++x) {
			for (int y = 0; y < 3; ++y) {
				for (int z = 0; z < 3; ++z) {
					gY += filterY[z][y][x] * Neighborhood[z][y][x].Get(3);
				}
			}
		}
		for (int x = 0; x < 3; ++x) {
			for (int y = 0; y < 3; ++y) {
				for (int z = 0; z < 3; ++z) {
					gX += filterX[z][y][x] * Neighborhood[z][y][x].Get(3);
				}
			}
		}

		//set after normalization
		OutputVoxel.Set(0, gX / 9.0f);
		OutputVoxel.Set(1, gY / 9.0f);
		OutputVoxel.Set(2, gZ / 9.0f);

		// END TODO 2
	};
};

EditorVisLabFilterExercise::EditorVisLabFilterExercise(Plugin & pluPlugin) : PluginInstance(pluPlugin),
		m_pTimerParameterChanged(NULL),
		m_bInputVolumeChanged(true),
		m_bOutputVolumeChanged(true),
		m_bParametersChanged(true),
		m_bComputationRunning(false),
		m_bAbortComputation(false),
		m_iKernelWidth(0),
		m_iKernelHeight(0),
		m_iKernelDepth(0)
{
	//Initialize parameter and observer
	m_modInputVolumeObserver.connect(this, &EditorVisLabFilterExercise::inputVolumeModified);
	m_modOutputVolumeObserver.connect(this, &EditorVisLabFilterExercise::outputVolumeModified);
	m_modParametersObserver.connect(this, &EditorVisLabFilterExercise::parameterModified);


	GetPlugin().GetProperty("Input Volume").require(Variant::TypeHandle()).addObserver(&m_modInputVolumeObserver);
	GetPlugin().GetProperty("Output Volume").require(Variant::TypeHandle()).addObserver(&m_modOutputVolumeObserver);

	GetPlugin().GetProperty("Smoothing Filter").require((Variant::TypeOption(),
													Variant("None"),
													Variant("Average"),
													Variant("Gaussian"),
													Variant("Bilateral")
													));

	GetPlugin().GetProperty("Gradient Filter").require((Variant::TypeOption(),
													Variant("None"),
													Variant("Backward Difference"),
													Variant("Central Difference"),
													Variant("Sobel")
													));

	GetPlugin().GetProperty("Enabled").addObserver(&m_modParametersObserver);
	GetPlugin().GetProperty("Auto-Recompute").require(Variant::TypeBoolean(true)).addObserver(&m_modParametersObserver);
	GetPlugin().GetProperty("Smoothing Filter").addObserver(&m_modParametersObserver);
	GetPlugin().GetProperty("Gradient Filter").addObserver(&m_modParametersObserver);

	GetPlugin().GetProperty("Gaussian Filter - Sigma").require(Variant(0.5f,0.5f,3.f));
	GetPlugin().GetProperty("Gaussian Filter - Sigma").addObserver(&m_modParametersObserver);

	GetPlugin().GetProperty("Bilateral Filter - Sigma Domain").require(Variant(0.5f,0.5f,3.f));
	GetPlugin().GetProperty("Bilateral Filter - Sigma Domain").addObserver(&m_modParametersObserver);

	GetPlugin().GetProperty("Bilateral Filter - Sigma Range").require(Variant(0.1f,0.1f,3.f));
	GetPlugin().GetProperty("Bilateral Filter - Sigma Range").addObserver(&m_modParametersObserver);

	//Timer for parameter
	//When it is changed wait a specific amount of time before starting calculation
	m_pTimerParameterChanged = new QTimer(this);
	m_pTimerParameterChanged->setInterval(500);		//Wait before start of calculation
	m_pTimerParameterChanged->setSingleShot(true);
	connect(m_pTimerParameterChanged,SIGNAL(timeout()),this,SLOT(triggerFilter()));

	m_pProgressWidget = new QFilterProgressWidget();
	m_pProgressWidget->setMaximumHeight(14);

	QHBoxLayout *pProgressWidgetLayout = new QHBoxLayout(m_pProgressWidget);
	pProgressWidgetLayout->setMargin(0);

	m_pProgressStatusLabel = new QLabel(m_pProgressWidget);
	m_pProgressStatusLabel->setText("");

	m_pProgressBar = new QProgressBar(m_pProgressWidget);
	// set some progress bar properties
	m_pProgressBar->setMinimum(0);
	m_pProgressBar->setMaximum(1000);
	m_pProgressBar->setValue(0);
	m_pProgressBar->setTextVisible(false);

	m_pProgressAbortButton = new QToolButton(m_pProgressWidget);
	m_pProgressAbortButton->setContentsMargins(0,0,0,0);
	m_pProgressAbortButton->setAutoRaise(true);
	m_pProgressAbortButton->setIcon(QApplication::style()->standardIcon(QStyle::SP_BrowserStop));
	m_pProgressAbortButton->setText("Abort");
	// connect some events to the button
	connect(m_pProgressAbortButton,SIGNAL(clicked()),this,SLOT(abortButtonClicked()));

	pProgressWidgetLayout->addWidget(m_pProgressStatusLabel);
	pProgressWidgetLayout->addWidget(m_pProgressBar);
	pProgressWidgetLayout->addWidget(m_pProgressAbortButton);

	m_pWatcher = new QFutureWatcher<void>();
	connect(m_pWatcher,SIGNAL(finished()),this,SLOT(computationFinished()));

	bool bEnabled = GetPlugin().GetProperty("Enabled");
	bool bRecompute = GetPlugin().GetProperty("Auto-Recompute");

	if (bEnabled && bRecompute)
		m_pTimerParameterChanged->start();
}

/**
 * @brief Destructor
 */
EditorVisLabFilterExercise::~EditorVisLabFilterExercise()
{
	if (m_pTimerParameterChanged)
	{
		delete m_pTimerParameterChanged;
		m_pTimerParameterChanged = NULL;
	}

	m_bAbortComputation = true;
	m_pWatcher->waitForFinished();

	if (m_pWatcher)
	{
		delete m_pWatcher;
		m_pWatcher = NULL;
	}

	if (m_pProgressWidget)
	{
		delete m_pProgressWidget;
		m_pProgressWidget = NULL;
	}
}

/**
 * @copydoc Editor::GetWidget()
 */
Widget * EditorVisLabFilterExercise::GetWidget()
{
	return NULL;
}

/**
 * @copydoc Editor::GetToolWidgets()
 */
std::list<Widget*> EditorVisLabFilterExercise::GetToolWidgets()
{
	std::list<Widget*> lisWidgets;
	return lisWidgets;
}

/**
 * @copydoc Editor::GetDockWidgets()
 */
std::list<Widget*> EditorVisLabFilterExercise::GetDockWidgets()
{
	return std::list<Widget*>();
}

/**
 * @copydoc Editor::GetMainActions()
 */
std::list<Action*> EditorVisLabFilterExercise::GetMainActions()
{
	return std::list<Action*>();
}

/**
 * @copydoc Editor::GetQuickActions()
 */
std::list<Action*> EditorVisLabFilterExercise::GetQuickActions()
{
	return std::list<Action*>();
}

/**
 * @copydoc Editor::onInitialized()
 */
void EditorVisLabFilterExercise::onInitialized()
{

}

/**
 * Event callback. This method is called when input volume has been changed
 **/
void EditorVisLabFilterExercise::inputVolumeModified(const Variant & varVariant, const Observable::Event & eveEvent)
{
	if (m_pTimerParameterChanged)
	{
		m_pTimerParameterChanged->stop();

		m_bAbortComputation = true;

		if (m_pWatcher)
			m_pWatcher->waitForFinished();

		m_bInputVolumeChanged = true;

		bool bEnabled = GetPlugin().GetProperty("Enabled");
		bool bRecompute = GetPlugin().GetProperty("Auto-Recompute");

		if (bEnabled && bRecompute)
			m_pTimerParameterChanged->start();
	}
}

/**
 * Event callback. This method is called when output volume has been changed
 **/
void EditorVisLabFilterExercise::outputVolumeModified(const Variant & varVariant, const Observable::Event & eveEvent)
{
	if (m_pTimerParameterChanged)
	{
		m_pTimerParameterChanged->stop();

		m_bAbortComputation = true;

		if (m_pWatcher)
			m_pWatcher->waitForFinished();

		m_bOutputVolumeChanged = true;

		bool bEnabled = GetPlugin().GetProperty("Enabled");
		bool bRecompute = GetPlugin().GetProperty("Auto-Recompute");

		if (bEnabled && bRecompute)
			m_pTimerParameterChanged->start();
	}
}

/** 
* Callback function for filter parameters
* When one of these parameters has been changed, start a timer. 
* Timer is necessary to prevent recalculation after each update of parameter. 
* The parameters are changed by slider control, and so you get a lot of updates when changing the slider.
* Recalculation is started when timer expires. 
*/ 
void EditorVisLabFilterExercise::parameterModified(const Variant & varVariant, const Observable::Event & eveEvent)
{
	if (m_pTimerParameterChanged)
	{
		m_pTimerParameterChanged->stop();

		m_bAbortComputation = true;

		if (m_pWatcher)
			m_pWatcher->waitForFinished();

		m_bParametersChanged = true;

		bool bEnabled = GetPlugin().GetProperty("Enabled");
		bool bRecompute = GetPlugin().GetProperty("Auto-Recompute");

		if (bEnabled && bRecompute)
			m_pTimerParameterChanged->start();
	}
}


void EditorVisLabFilterExercise::abortButtonClicked()
{
	if (m_pTimerParameterChanged)
		m_pTimerParameterChanged->stop();

	m_bAbortComputation = true;

	if (m_pWatcher)
		m_pWatcher->waitForFinished();


	computationProgress();
}

void EditorVisLabFilterExercise::triggerFilter()
{
	startComputation();				
}

void EditorVisLabFilterExercise::startComputation()
{
	if (!m_bParametersChanged && !m_bInputVolumeChanged && !m_bOutputVolumeChanged)
		return;

	if (m_bComputationRunning)
		return;

	Handle hanInput = GetPlugin().GetProperty("Input Volume");
	VolumeResource *pInput  = hanInput.GetResource<VolumeResource>();

	if (!pInput)
		return;

	Handle hanOutput = GetPlugin().GetProperty("Output Volume");
	VolumeResource *pOutput  = hanOutput.GetResource<VolumeResource>();

	if (!pOutput)
		return;

	m_bParametersChanged = false;
	m_bAbortComputation = false;
	m_bComputationRunning = true;

	m_strSmoothingFilter = GetPlugin().GetProperty("Smoothing Filter");
	m_strGradientFilter = GetPlugin().GetProperty("Gradient Filter");
	m_fGaussianSigma = GetPlugin().GetProperty("Gaussian Filter - Sigma");
	m_fBilateralSigmaD = GetPlugin().GetProperty("Bilateral Filter - Sigma Domain");
	m_fBilateralSigmaR = GetPlugin().GetProperty("Bilateral Filter - Sigma Range");

	m_fProgress = 0.0f;
	m_strMessage = "Starting";
	computationProgress();

	GetPlugin().addStatusWidget(m_pProgressWidget);

	m_volFilter.clear();
	m_volFilter.copyProperties(*pInput);
	m_volFilter.assign(*pInput);

	QFuture<void> future = QtConcurrent::run(this,&EditorVisLabFilterExercise::performComputation);
	m_pWatcher->setFuture(future);
}


void EditorVisLabFilterExercise::performComputation()
{
	if(m_strSmoothingFilter == "Average")
	{
		m_strMessage = "Smoothing (Average)";
		computeAverageKernel(3,3,3);

		// filtering function with iterator
		VOLUME_CALL(computeConvolution, m_volFilter);
	}
	else if (m_strSmoothingFilter == "Gaussian")
	{
		m_strMessage = "Smoothing (Gaussian)";
		computeGaussianKernel(m_fGaussianSigma);

		// filtering function with iterator
		VOLUME_CALL(computeConvolution, m_volFilter);
	}
	else if (m_strSmoothingFilter == "Bilateral")
	{
		m_strMessage = "Smoothing (Bilateral)";
		computeGaussianKernel(m_fBilateralSigmaD);

		// filtering function with iterator - has to be implemented
		VOLUME_CALL(computeBilateral, m_volFilter,m_fBilateralSigmaR);
	}

	if (!m_bAbortComputation)
	{
		if(m_strGradientFilter == "Backward Difference")
		{
			m_strMessage = "Computing Gradients (Backward Differences)";
			computeGradients< BackwardDifferenceGradientFilter >(); 
		}
		else if (m_strGradientFilter == "Central Difference")
		{
			m_strMessage = "Computing Gradients (Central Differences)";
			computeGradients< CentralDifferenceGradientFilter >(); 
		}
		else if (m_strGradientFilter == "Sobel")
		{
			m_strMessage = "Computing Gradients (Sobel)";
			computeGradients< SobelGradientFilter >(); 
		}
	}
}

void EditorVisLabFilterExercise::computationProgress()
{
	if (m_pProgressWidget)
	{
		m_pProgressStatusLabel->setText(m_strMessage);
		m_pProgressBar->setValue(m_fProgress*m_pProgressBar->maximum());
	}
}

void EditorVisLabFilterExercise::computationFinished()
{
	m_strMessage = "Finishing";
	computationProgress();

	Handle hanOutput = GetPlugin().GetProperty("Output Volume");
	VolumeResource *pOutput  = hanOutput.GetResource<VolumeResource>();

	if (pOutput)
	{
		if (!m_bAbortComputation)
			VOLUME_CALL(copyResult,*pOutput);
	}

	GetPlugin().removeStatusWidget(m_pProgressWidget);

	m_bComputationRunning = false;
	m_bAbortComputation = false;

}

template<template <class,unsigned int> class FunctionType>
void EditorVisLabFilterExercise::computeGradients()
{
	FunctionType<float,4> filterFunction;
	NeighborhoodTransformer<FunctionType,Volume,float,4> transformer(filterFunction, m_volFilter, m_volFilter);

	while (!transformer.IsAtEnd())
	{
		transformer.transform();
		m_fProgress = transformer.GetProgress();
		QMetaObject::invokeMethod(this,"computationProgress",Qt::QueuedConnection);

		if (m_bAbortComputation)
		{
			transformer.cancel();
			return;
		}
	}
}

// TODO 3: Implement 3D Convolution (convolution is only in 2D slices right now)
// http://en.wikipedia.org/wiki/Convolution#Discrete_convolution
// Matlab in 2D convolution
// http://www.mathworks.de/help/techdoc/ref/conv2.html

template <class ElementType,unsigned int N>
void EditorVisLabFilterExercise::computeConvolution(const Volume<ElementType,N> & volVolume)
{
	Volume<ElementType,N> volSrc = volVolume;

	Volume<ElementType,N>::BlockIterator blockIter(volSrc);
	Volume<float,4>::LinearManipulator manip(m_volFilter);
	
	int size_rows = volSrc.GetWidth();
	int size_columns = volSrc.GetHeight();
	int size_slices = volSrc.GetDepth();

	unsigned int uVoxelCount = volSrc.GetWidth()*volSrc.GetHeight()*volSrc.GetDepth();
	unsigned int uVoxelIndex = 0;

	// maybe should be member of class
	int kernel_size_rows = m_iKernelWidth;
	int kernel_size_columns = m_iKernelHeight;
	int kernel_size_slices = m_iKernelDepth;

	///Do convolution in X, Y and Z direction
	while (!blockIter.IsAtEnd())
	{
		unsigned int uBlockX = blockIter.GetPositionX()*BLOCKDIMENSION;
		unsigned int uBlockY = blockIter.GetPositionY()*BLOCKDIMENSION;
		unsigned int uBlockZ = blockIter.GetPositionZ()*BLOCKDIMENSION;

		Volume<ElementType,N>::BlockLinearIterator voxelIter(blockIter);

		while (!voxelIter.IsAtEnd())
		{
			unsigned int uVoxelX = uBlockX+voxelIter.GetPositionX();
			unsigned int uVoxelY = uBlockY+voxelIter.GetPositionY();
			unsigned int uVoxelZ = uBlockZ+voxelIter.GetPositionZ();

			float fConvoluted = 0.0f;

			int offset_columns = (-1)*((kernel_size_columns - 1)/2);

			for (int j = 0; j < kernel_size_columns; j++)
			{
				int offset_rows = (-1)*((kernel_size_rows - 1)/2);

				for (int i = 0; i < kernel_size_rows; i++)
				{
					Voxel<float,1> voxValue = voxelIter.GetNeighbor(offset_rows,offset_columns,0);
					float signal_value = voxValue.Get(0);

					int linear_index = j*(kernel_size_rows) + i;
					fConvoluted = fConvoluted + signal_value * m_vecConvolutionKernel[linear_index];
					offset_rows = offset_rows + 1;
				}

				offset_columns = offset_columns + 1;
			}

			manip.SetPosition(uVoxelX,uVoxelY,uVoxelZ);

			Voxel<float,4> voxel = manip.Get();
			voxel.Set(3,fabsf(fConvoluted));
			manip.Set(voxel);

			++uVoxelIndex;
			++voxelIter;
		}

		// set new normalized progress
		m_fProgress = float(uVoxelIndex)/float(uVoxelCount);
		QMetaObject::invokeMethod(this,"computationProgress",Qt::QueuedConnection);

		if (m_bAbortComputation)
			break;

		++blockIter;
	}
}
// END TODO 3

// TODO BONUS: Implement 3D Bilateral filter
//
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.4.634&rep=rep1&type=pdf
//
// http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/MANDUCHI1/Bilateral_Filtering.html
// 
template <class ElementType,unsigned int N>
void EditorVisLabFilterExercise::computeBilateral(const Volume<ElementType,N> & volVolume, 	float fSigmaR)
{
}
// END TODO BONUS


template <class ElementType,unsigned int N>
void EditorVisLabFilterExercise::copyResult( Volume<ElementType,N> & volResult )
{
	m_modOutputVolumeObserver.SetEnabled(false);
	volResult.assign(m_volFilter);
	volResult.copyProperties(m_volFilter);
	m_modOutputVolumeObserver.SetEnabled(true);
}

// http://en.wikipedia.org/wiki/Moving_average
void EditorVisLabFilterExercise::computeAverageKernel(int iWidth, int iHeight, int iDepth)
{
	// alloc memory for 3D kernel
	m_vecConvolutionKernel.assign(iWidth*iHeight*iDepth,0.0f);

	int x_offset = (-1)*(iWidth-1)/2;
	int y_offset = (-1)*(iHeight-1)/2;
	int z_offset = (-1)*(iDepth-1)/2;

	float norm_sum = 0.0f;

	for(int k = 0; k < iDepth; k++)
	{
		y_offset = (-1)*(iHeight-1)/2;
		for(int j = 0; j < iHeight; j++)
		{
			x_offset = (-1)*(iWidth-1)/2;
			for(int i = 0; i < iWidth; i++)
			{
				int linear_index = k*(iWidth*iHeight) + j*(iWidth) + i;

				// 3d average filter
				m_vecConvolutionKernel[linear_index] = 1;
				// sum up noramlization denominator
				norm_sum = norm_sum + m_vecConvolutionKernel[linear_index];

				x_offset = x_offset + 1;
			}
			y_offset = y_offset + 1;
		}
		z_offset = z_offset + 1;
	}

	// normalize to one
	for(int i = 0; i < m_vecConvolutionKernel.size(); i++)
	{
		m_vecConvolutionKernel[i] = m_vecConvolutionKernel[i] / norm_sum;
	}

	m_iKernelWidth = iWidth;
	m_iKernelHeight = iHeight;
	m_iKernelDepth = iDepth;
}

void EditorVisLabFilterExercise::computeGaussianKernel( float fSigma )
{
	// TODO 4: Compute the correct minimum size of the Gaussian smoothing filter kernel
	// gaussian kernel is (6*sigma - 1) in one dimension
	// http://en.wikipedia.org/wiki/Gaussian_filter
	int iKernelWidth = 3;
	// END TODO 4

	// we make it alway at least 3x3x3 size
	if( iKernelWidth < 3)
		iKernelWidth = 3;

	computeGaussianKernel(fSigma,iKernelWidth,iKernelWidth,iKernelWidth);
}

void EditorVisLabFilterExercise::computeGaussianKernel(float fSigma, int iWidth, int iHeight, int iDepth )
{
	m_vecConvolutionKernel.assign(iWidth*iHeight*iDepth,0.0f);

	int x_offset = (-1)*(iWidth-1)/2;
	int y_offset = (-1)*(iHeight-1)/2;
	int z_offset = (-1)*(iDepth-1)/2;

	float norm_sum = 0.0f;

	for(int k = 0; k < iDepth; k++)
	{
		y_offset = (-1)*(iWidth-1)/2;
		for(int j = 0; j < iWidth; j++)
		{
			x_offset = (-1)*(iWidth-1)/2;
			for(int i = 0; i < iWidth; i++)
			{
				int linear_index = k*(iWidth*iWidth) + j*(iWidth) + i;

				// TODO 5: Fill a Gaussian kernel with real values
				// http://en.wikipedia.org/wiki/Gaussian_filter
				// 3d gaussian function
				m_vecConvolutionKernel[linear_index] = 1;
				// END TODO 5

				norm_sum = norm_sum + m_vecConvolutionKernel[linear_index];

				x_offset = x_offset + 1;
			}
			y_offset = y_offset + 1;
		}
		z_offset = z_offset + 1;
	}

	// normalize to one
	for(int i = 0; i < m_vecConvolutionKernel.size(); i++)
	{
		m_vecConvolutionKernel[i] = m_vecConvolutionKernel[i] / norm_sum;
	}

	m_iKernelWidth = iWidth;
	m_iKernelHeight = iHeight;
	m_iKernelDepth = iDepth;
}


