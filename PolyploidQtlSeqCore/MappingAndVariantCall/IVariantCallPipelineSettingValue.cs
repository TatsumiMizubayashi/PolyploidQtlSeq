using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.VariantCall;

namespace PolyploidQtlSeqCore.MappingAndVariantCall
{
    /// <summary>
    /// 変異検出パイプライン設定値インターフェース
    /// </summary>
    public interface IVariantCallPipelineSettingValue : IMappingSettingValue, IMappingSampleSettingValue, IVariantCallScenarioSettingValue
    {
    }
}
