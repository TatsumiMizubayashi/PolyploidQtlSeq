using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.VariantCall;

namespace PolyploidQtlSeqCore.MappingAndVariantCall
{
    /// <summary>
    /// 変異検出パイプライン設定
    /// </summary>
    internal class VariantCallPipelineSettings
    {
        /// <summary>
        /// 変異検出パイプライン設定インスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">設定値</param>
        public VariantCallPipelineSettings(IVariantCallPipelineSettingValue settingValue)
        {
            MappingSettings = new MappingSettings(settingValue);
            MappingSampleSettings = new MappingSampleSettings(settingValue);
            VariantCallScenarioSettings = new VariantCallScenarioSettings(settingValue);
        }

        /// <summary>
        /// Mapping設定を取得する。
        /// </summary>
        public MappingSettings MappingSettings { get; }

        /// <summary>
        /// Mappingサンプル設定を取得する。
        /// </summary>
        public MappingSampleSettings MappingSampleSettings { get; }

        /// <summary>
        /// 変異検出シナリオ設定を取得する。
        /// </summary>
        public VariantCallScenarioSettings VariantCallScenarioSettings { get; }
    }
}
