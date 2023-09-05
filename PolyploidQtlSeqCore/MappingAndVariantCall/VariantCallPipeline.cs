using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.VariantCall;

namespace PolyploidQtlSeqCore.MappingAndVariantCall
{
    /// <summary>
    /// 変異検出パイプライン
    /// </summary>
    internal class VariantCallPipeline
    {
        private readonly VariantCallPipelineSettings _settings;

        /// <summary>
        /// 変異検出パイプラインインスタンスを作成する。
        /// </summary>
        /// <param name="settings">設定</param>
        public VariantCallPipeline(VariantCallPipelineSettings settings)
        {
            _settings = settings;
        }

        /// <summary>
        /// パイプラインを実行する。
        /// </summary>
        /// <returns>VCFファイル</returns>
        public async ValueTask<VcfFile> RunAsync()
        {
            var allSampleBamFiles = await MappingAsync();
            return await VariantCallAsync(allSampleBamFiles);
        }

        /// <summary>
        /// 全サンプルのMappingを行う。
        /// </summary>
        /// <returns>全サンプルのBAMファイル</returns>
        private async ValueTask<AllSampleBamFiles> MappingAsync()
        {
            var mappingScenario = new AllSampleMappingScenario(_settings.MappingSettings);
            return await mappingScenario.MappingAsync(_settings.MappingSampleSettings);
        }

        /// <summary>
        /// 変異検出を行う。
        /// </summary>
        /// <param name="allBams">全サンプルのBAMファイル</param>
        /// <returns>VCFファイル</returns>
        private async ValueTask<VcfFile> VariantCallAsync(AllSampleBamFiles allBams)
        {
            var variantCallScenario = new VariantCallScenario(_settings.VariantCallScenarioSettings);
            return await variantCallScenario.CallAsync(allBams);
        }
    }
}
