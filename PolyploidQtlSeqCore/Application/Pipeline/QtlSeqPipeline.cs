using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.Application.QtlAnalysis;
using PolyploidQtlSeqCore.IO;
using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.Mapping;
using PolyploidQtlSeqCore.QtlAnalysis.VariantCall;
using PolyploidQtlSeqCore.VariantCall;

namespace PolyploidQtlSeqCore.Application.Pipeline
{
    /// <summary>
    /// QTL-Seqパイプライン
    /// </summary>
    public class QtlSeqPipeline
    {
        private readonly QtlSeqCommandOption _option;
        private readonly OutputDirectory _outputDir;
        private readonly ThreadNumber _threadNumber;

        /// <summary>
        /// QTL-Seqパイプラインインスタンスを作成する。
        /// </summary>
        /// <param name="optionValues">オプションの値</param>
        /// <param name="options">コマンドオプション</param>
        public QtlSeqPipeline(IQtlSeqCommandOptions optionValues, IReadOnlyCollection<CommandOption> options)
        {
            _option = new QtlSeqCommandOption(optionValues, options);
            _outputDir = _option.QtlAnalysisScenarioOptions.OutputDir;
            _threadNumber = _option.QtlAnalysisScenarioOptions.ThreadNumber;
        }

        /// <summary>
        /// QTL-Seq解析を実行する。
        /// </summary>
        /// <returns>終了コード</returns>
        public async ValueTask<int> RunAsync()
        {
            var code = 0;

            try
            {
                CommandLog.Clear();
                _outputDir.Create();

                var allSampleBamFiles = await MappingAsync();
                var vcfFile = await VariantCallAsync(allSampleBamFiles);
                code = QtlAnalysis(vcfFile);
            }
            catch (Exception ex)
            {
                code = 1;
                Console.Error.WriteLine(ex.Message);
                Console.Error.WriteLine(ex.StackTrace);
            }
            finally
            {
                var commandLogFilePath = _outputDir.CreateFilePath("Command Log.txt");
                CommandLog.Save(commandLogFilePath);

                var paramsFilePath = _outputDir.CreateFilePath("qtlSeq.parameter.txt");
                _option.SaveParameterFile(paramsFilePath);
            }

            return code;
        }

        private async ValueTask<AllSampleBamFiles> MappingAsync()
        {
            var allSampleMappingScenario = new AllSampleMappingScenario(_option.ReferenceSequence, _threadNumber);
            return await allSampleMappingScenario.MappingAsync(_option.MappingOption).AsTask();
        }

        private ValueTask<VcfFile> VariantCallAsync(AllSampleBamFiles allSampleBamFiles)
        {
            var analysisChrs = allSampleBamFiles.GetAnalysisChrs(_option.AnalysisChrOption);

            var variantCallScenario = new VariantCallScenario(
                _option.BcfToolsVariantCallOption,
                _option.SnpEffOption,
                _option.ReferenceSequence,
                _threadNumber);

            return variantCallScenario.CallAsync(allSampleBamFiles, _outputDir, analysisChrs);
        }

        private int QtlAnalysis(VcfFile vcfFile)
        {
            var dummyParameterDictionary = new Dictionary<string, string>();
            var dummyUserOptionDictionary = new Dictionary<string, bool>();
            var inputVcf = new InputVcf(vcfFile.Path, dummyParameterDictionary, dummyUserOptionDictionary);

            var qtlAnalysisScenario = new QtlAnalysisScenario(_option.QtlAnalysisScenarioOptions);
            return qtlAnalysisScenario.Run(inputVcf);
        }
    }
}
