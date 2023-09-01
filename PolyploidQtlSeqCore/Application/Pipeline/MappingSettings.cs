using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.Share;

namespace PolyploidQtlSeqCore.Application.Pipeline
{
    /// <summary>
    /// Mapping設定
    /// </summary>
    internal class MappingSettings
    {
        /// <summary>
        /// Mapping設定インスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">設定値</param>
        public MappingSettings(IMappingSettingValue settingValue)
        {
            ReferenceSequence = new ReferenceSequence(settingValue.ReferenceSequence);
            ThreadNumber = new ThreadNumber(settingValue.ThreadNumber);
        }

        /// <summary>
        /// リファレンスシークエンス
        /// </summary>
        public ReferenceSequence ReferenceSequence { get; }

        /// <summary>
        /// スレッド数を取得する。
        /// </summary>
        public ThreadNumber ThreadNumber { get; }
    }
}
