using PolyploidQtlSeqCore.IO;
using PolyploidQtlSeqCore.MappingAndVariantCall;
using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.VariantCall;

namespace PolyploidQtlSeqCore.Application.Pipeline
{
    /// <summary>
    /// QTL-Seqパイプライン
    /// </summary>
    public class QtlSeqPipeline
    {
        private readonly VariantCallPipelineSettings _variantCallPipelineSettings;
        private readonly QtlAnalysisScenarioSettings _qtlAnalysisScenarioSettings;

        /// <summary>
        /// QTL-Seqパイプラインインスタンスを作成する。
        /// </summary>
        /// <param name="optionValues">オプションの値</param>
        [Obsolete("削除予定")]
        public QtlSeqPipeline(IQtlSeqPipelineSettingValue optionValues)
        {
            _variantCallPipelineSettings = new VariantCallPipelineSettings(optionValues);
            _qtlAnalysisScenarioSettings = new QtlAnalysisScenarioSettings(optionValues);
        }

        /// <summary>
        /// QTL-Seqパイプラインインスタンスを作成する。
        /// </summary>
        /// <param name="variantCallPipelineSettingValue">変異検出パイプライン設定値</param>
        /// <param name="qtlAnalysisScenarioSettingValue">QTL解析設定値</param>
        public QtlSeqPipeline(IVariantCallPipelineSettingValue variantCallPipelineSettingValue,
            IQtlAnalysisScenarioSettingValue qtlAnalysisScenarioSettingValue)
        {
            _variantCallPipelineSettings = new VariantCallPipelineSettings(variantCallPipelineSettingValue);
            _qtlAnalysisScenarioSettings = new QtlAnalysisScenarioSettings(qtlAnalysisScenarioSettingValue);
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

                var vcfFile = await VariantCallAsync();
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
                var commandLogFilePath = _qtlAnalysisScenarioSettings.OutputDir.CreateFilePath("Command Log.txt");
                CommandLog.Save(commandLogFilePath);
            }

            return code;
        }


        private async ValueTask<VcfFile> VariantCallAsync()
        {
            var pipeline = new VariantCallPipeline(_variantCallPipelineSettings);
            return await pipeline.RunAsync();
        }

        private int QtlAnalysis(VcfFile vcfFile)
        {
            var inputVcf = new InputVcf(vcfFile.Path);

            var qtlAnalysisScenario = new QtlAnalysisScenario(_qtlAnalysisScenarioSettings);
            return qtlAnalysisScenario.Run(inputVcf);
        }
    }
}
