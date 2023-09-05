using PolyploidQtlSeqCore.QtlAnalysis.Chr;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// 変異検出シナリオ設定値インターフェース
    /// </summary>
    public interface IVariantCallScenarioSettingValue : IAnalysisChrSettingValue, IBcftoolsVariantCallSettingValue, ISnpEffSettingValue
    {
    }
}
